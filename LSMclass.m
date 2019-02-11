classdef LSMclass <handle
    %LSMObject constructor  
    %requires a list of images (full file paths), file type, and image resize rate 
    % LSMObject(listOfFiles,fileType,resizeRate)
    %   Detailed explanation goes here
    
    properties (Access=protected)
        imds
        fileType='png'
        phi
        resizeRate
        initiArgs
        lsmArgs
        contours
    end
    
    methods 

        function this = LSMObject(imds)
            this.imds= imds; %dir(imds.inputPath);
            %this.inFolder=this.imds.folder;
            this.resizeRate=.5;
            this.setLSMOptions()
            
            % Set file type
            %this.fileType=fileType;
            
            % Initialize Contours
            this.contours=cell(numel(this.imds),...
                1+ceil(this.lsmArgs.outerIter*this.lsmArgs.innerIter/10));
        end
        
        function setResize(this, rate)
            this.resizeRate=rate;
        end
        
        function  setLSMOptions(this, varargin)
            try
                parser = createParserLSM();
                parser.parse(varargin{:});
                [this.initiArgs, this.lsmArgs] = convertToCanonicalFormLSM(parser);
                this.initiArgs.resize=this.resizeRate;
                this.contours=cell(numel(this.imds),...
                1+ceil(this.lsmArgs.outerIter*this.lsmArgs.innerIter/10));
            catch e
                % Reduce the stack trace of the error message by throwing as caller
                throwAsCaller(e)
            end
        end
        
        function out = getLSMoptions(this)
            out.lambda = this.lsmArgs.lambda;
            out.alpha = this.lsmArgs.alpha;
            out.mu = this.lsmArgs.mu;
            out.innerIter = this.lsmArgs.innerIter;
            out.outerIter = this.lsmArgs.outerIter;
            out.timestep = this.lsmArgs.timestep;
            out.narrowBand = this.lsmArgs.narrowBand;
            out.snapshots = this.lsmArgs.snapshots;
        end
        
        
        function cs=getContours(this,imgIdx,step)
            cs= this.contours{imgIdx,step};
        end
        
        function last=lastContour(this)
            last= size(this.contours,2);
        end
        
        function [fname,ftype]=fileInfo(this,k)
              [fname,ftype]= getFileName(this.imds{k});
        end
        
        function showContours(this, k,s)
                [img, ~, ftype]=readImage(this.imds{k});
                ce=this.contours{k,s}; %returns the contours for image k
                showImage(img, ftype)  %show the image   
                hold on
                for j=1:numel(ce) %draw the contour
                    plot(ce{j}(1,:), ce{j}(2,:), 'b');
                end       
        end
        
        function [fname,ftype]=disp(this,k)
            [fname,ftype]=getFileName(this.imds{k});
            this.showContours(k,this.lastContour);
        end
        
        function runLSM(this)
            for k=1:numel(this.imds)
                h=figure(k);
                %read the image
                [img,~,~]=readImage(this.imds{k});
                img= imresize(img, this.resizeRate);
                g=indicateEdge(img);
                %set initial lsf and get the 
                this.phi=initializeLSM(img, this.initiArgs);
                this.contours{k,1}= getLSF(this.phi, 1/this.resizeRate);
                dispc(h, img, this.phi,this.imds{k}, 'initial');
                % LSM
                for n=1:this.lsmArgs.outerIter
                    this.phi=lsmReg(this.phi, g, this.lsmArgs);
                    if rem(n,this.lsmArgs.snapshots)==0
                        dispc(h, img, this.phi, this.imds{k},  num2str(this.lsmArgs.innerIter*n));
                    end
                    m=this.lsmArgs.innerIter*n;
                    if rem(m,10)==0
                        this.contours{k,m/10+1}= getLSF(this.phi, 1/this.resizeRate);
                    end
                end
                 this.contours{k,end}= getLSF(this.phi, 1/this.resizeRate);
                 dispc(h, img, this.phi, this.imds{k},  [num2str(this.lsmArgs.innerIter*n) '_final']);   
            end
        end
    end
end
%% initialization
 function phi=initializeLSM(img, initiArgs)
    c0=initiArgs.c0; 
    x=(initiArgs.x).*(initiArgs.resize); dx=initiArgs.dx; 
    y=(initiArgs.y).*(initiArgs.resize); dy=initiArgs.dy;
    initialLSF = c0*ones(size(img));
    for i=1:numel(x)
    if x(i)>0
    initialLSF(:, x(i):x(i)+dx)=-c0; 
    end
    if y(i)>0
    initialLSF(y(i):y(i)+dy, :)=-c0;
    end
    end
    phi=initialLSF;            
 end
%% LSF and Edge         
 function lsc=getLSF(phi, rs)
    phi=imresize(phi,rs);
    c = contourc(phi, [0,0]);
    x1=1;
    k=0;
    lsc=cell(1,100);
    while x1<size(c,2)
        k=k+1;
        x2=c(2,x1)+x1+1;
        lsc{k}=[c(1,x1+1:x2-1) ; c(2,x1+1:x2-1)];
        x1=x2;
    end
    if k>0 && k<100
        lsc=lsc(1:k);
    end
 end
         
 function  g=indicateEdge(image)
    G=fspecial('average',5); % Average kernel
    img_smooth=conv2(image,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.
 end
 %% show progress
 
 function dispc(h,image, phi, filepath,  strtitle) %display image and contours 
    figure(h);
    [fname,ftype]=getFileName(filepath);
    showImage(image,ftype);
    hold on;
    contour(phi, [0,0], 'r');
    %set (h2, 'LineWidth', 2);  
    title([fname , ': step ' , strtitle ], 'Interpreter', 'none');
 end
 
 function [fname,ftype]=getFileName(file)
        [a,~]=strtok(wrev(file), filesep);
        [fname,ftype]=strtok(wrev(a), '.');
 end
 %% parser
         
function p = createParserLSM()

p = inputParser;

%LSM Args
defaultLambda = 5;
defaultAlpah = -3;
defaultMu = 0.1;
defaultInnerIter = 2;
defaultOuterIter = 400;
defaultTimestep = 2;
defaultViewSnaps = 100;
defaultNarrowBand = false;

p.addParameter('lambda', defaultLambda, @numeric);
p.addParameter('alpha', defaultAlpah, @numeric);
p.addParameter('mu', defaultMu, @numeric);
p.addParameter('innerIter', defaultInnerIter, @isPositiveInteger);
p.addParameter('outerIter', defaultOuterIter, @isPositiveInteger);
p.addParameter('timestep', defaultTimestep, @isPositiveInteger);
p.addParameter('snapshots', defaultViewSnaps, @isPositiveInteger);
p.addParameter('narrowBand', defaultNarrowBand, @islogical);

% initial Args
defaultc0 = 2;
defaultx =0;
defaultdx = 5;
defaulty = 120;
defaultdy = 5;


p.addParameter('c0', defaultc0, @isnumeric);
p.addParameter('x', defaultx, @isnumeric);
p.addParameter('dx', defaultdx, @isnumeric);
p.addParameter('y', defaulty, @isnumeric);
p.addParameter('dy', defaultdy, @isnumeric);
end

function [initiArgs , lsmArgs] = convertToCanonicalFormLSM(parser)
results = parser.Results;
lsmArgs = struct;
lsmArgs.lambda = results.lambda;
lsmArgs.alpha = results.alpha;
lsmArgs.mu = results.mu;
lsmArgs.innerIter = results.innerIter;
lsmArgs.outerIter = results.outerIter;
lsmArgs.timestep = results.timestep;
lsmArgs.snapshots = results.snapshots;
lsmArgs.narrowBand = results.narrowBand;

initiArgs=struct;
initiArgs.c0 = results.c0;
initiArgs.x= results.x;
initiArgs.dx = results.dx;
initiArgs.y = results.y;
initiArgs.dy = results.dy;
end


function tf = isPositiveInteger(x)
isPositive = x>0;
isInteger = isreal(x) && isnumeric(x) && all(mod(x,1)==0);
tf = isPositive && isInteger;
end

function tf = numeric(x)
tf= isreal(x) && isnumeric(x) && isscalar(x);
end

%% lsmReg
function phi = lsmReg(phi0, g,lsmArgs )
%This Matlab code implements an edge-based active contour model without
%narrow band
%
%  Input:
%      phi0: level set function to be updated by level set evolution
%      g: edge indicator function
% .    lambda:  length termcoefficient
%      mu: distance regularization term coefficient
%      alfa:   area term coefficient
%      timestep: time step
%      
%      
%  Output:
%      phi: updated level set function after level set evolution
 
    lambda=lsmArgs.lambda;
    mu=lsmArgs.mu;
    alpha=lsmArgs.alpha;
    timestep=lsmArgs.timestep;
    iter=lsmArgs.innerIter;
    nb=lsmArgs.narrowBand;
    epsiln=1e-10;
    sigma=1.5; %width of Dirac Delta function
    phi=phi0;
    [vx, vy]=gradient(g);
    for k=1:iter %% iteration_inner
        phi=nbc(phi);
        phi1=phi;%nb
        [phi_x,phi_y]=gradient(phi);
        if nb
        phi=phi.*(phi_x+phi_y~=0); %nb
        end
        s=sqrt(phi_x.^2 + phi_y.^2);   
        Nx=phi_x./(s+epsiln); 
        Ny=phi_y./(s+epsiln);
        curvature=divergence(Nx,Ny);

        diracPhi=ddirac(phi,sigma);
    % terms
        areaTerm=diracPhi.*g; 
        edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
        regTerm=distReg(phi); 
        phi=phi1 + timestep*(mu*regTerm + lambda*edgeTerm + alpha*areaTerm);
    end
end
 
function f = ddirac(x, sigma)
    f=(1/2/sigma)*(1+cos(pi*x/sigma));
    b = (x<=sigma) & (x>=-sigma);
    f = f.*b;
end

function g = nbc(f)
    g = f;
    g(1,:)=f(3,:);
    g(end,:)=f(end-2,:);
    g(:,end)=g(:,end-2);
    g(:,1)=g(:,3);
end
