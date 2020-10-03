%===============================%
% Manhattan Lattice             %
% Addison Hartman               %
% Weld Lab at UCSB              %
% Last edited in October 2020   %
%===============================%

%% README

% To run a simulation, only edit values in "File outputs" and "Simulation
% Parameters". Know what you're doing before you edit anything else.

%% File Outputs

makeVfreeContourPlot = false;
makeVfreeSurfacePlot = false;
makeInitialStatePlot = false;
% to perfect appearance of figures for outside use, save the pop-up as a 
% .fig and use the "property inpector" or "export setup" tools to refine 
% axes labels, tick marks, etc.

%% Simulation Parameters

% space
Nx      = 2^11;     % number of points (make power of 2)
Ny      = 2^11;
dx      = 0.08;     % position mesh spacing in 1/kL
dy      = 0.08;

% saves (R for ramp, F for free evolution)
Rsaves      = 40;   % number of times to save the wave function
Fsaves      = 100;

% time (R for ramp, F for free evolution)
dtR_ms      = 0.01;	% time spacing in ms, make <0.00758 for dt<0.1 1/omegar
dtF_ms      = 0.01;
tRamp_ms    = 80;
tFree_ms    = 100;

%% Fundamental Constants

amu     = 1.66E-27;     % atomic mass unit in kg
h       = 6.626E-34;	% Planck's constant in J*s
hbar	= h/(2*pi);     % reduced Planck's constant in J*s

%% Experiment Parameters

m       = 84*amu;       % mass of strontium in kg
lambda  = 1064E-9;      % wavelength in m

omegax  = 270;          % trap frequencies from Rajagopal's thesis
omegay  = 430;

V0      = 5;            % lattice depth in ER
FX      = 0;            % x-direction force
FY      = -1*10^(-3);       % y-direction force

%% Derived Quantities

d       = lambda/2;             % lattice spacing in m
kL      = 2*pi/lambda;          % wave number in 1/m
Er      = (hbar*kL)^2/(2*m);    % recoil energy in J
omegar  = Er/hbar;              % recoil frequency in Hz

% natural units of the system
% length: 1/kL
% energy: Er
% time:   1/omegar

%% Position-Momentum Mesh

xVals   = dx*(-Nx/2:1:Nx/2-1);      yVals   = dy*(-Ny/2:1:Ny/2-1);  % position mesh points

dkx     = 2*pi/dx/(Nx-1);           dky     = 2*pi/dy/(Ny-1);       % k mesh spacing in kL
kxVals  = dkx*(-Nx/2:1:Nx/2-1);     kyVals  = dky*(-Ny/2:1:Ny/2-1); % momentum mesh points

[X,Y]   = meshgrid(xVals,yVals);    % position mesh used later
[Kx,Ky] = meshgrid(kxVals,kyVals);  % momentum mesh used later

% note fftshift() is used on the momentum operator (Uk -> UkShf) in place
% of the momentum matrix itself for efficiency

xMax    = Nx/2*dx;      % used in plots, units of 1/kL
yMax    = Ny/2*dx;
kxMax   = Nx/2*dkx;     % used in plots, units of
kyMax   = Nx/2*dky;

%% Time Discretization

% varR indicates during ramp-up of lattice
% varF indicates during free evolution

dtR         = dtR_ms/1000*omegar;   % time spacing in 1/omegar
dtF         = dtF_ms/1000*omegar;
tRamp       = tRamp_ms/1000*omegar;	% ramp time in 1/omegar
tFree       = tFree_ms/1000*omegar;	% final time in 1/omegar
tRsteps     = round(tRamp/dtR);     % number of time steps
tFsteps     = round(tFree/dtF);
tRamp       = dtR*tRsteps;          % re-define end time
tFree       = dtF*tFsteps;

Rsavestep   = tRsteps/(Rsaves-1);   % how often to save, accounting for initial wavefctn
Fsavestep   = tFsteps/Fsaves;
Rsave       = 1;                    % save counter
Fsave       = 0;                    % save counter

%% Potentials

Vlat    = (V0/2)*(cos(2.*X)+cos(2.*Y)+4*cos(X).*cos(Y));        % lattice potential
Vgrav   = -FX.*X-FY.*Y;                                         % force from lattice tilt
Vodt    = 1/2*m/(Er*kL^2).*((omegax^2.*X.^2)+(omegay^2.*Y.^2)); % optical dipole trap potential

Vt      = @(t) t/tRamp;                     % ramp-up factor
Vramp   = @(t) Vodt+Vt(t).*Vlat;            % total potential during ramp period
Vfree   = Vlat+Vgrav;                       % total potential during free evolution

% plots of Vfree (set to true in "File Outputs" to create)
if makeVfreeSurfacePlot
    makePotentialSurfacePlot(X,Y,Vfree,xMax,yMax) %#ok<UNRCH>
end
if makeVfreeContourPlot
    makePotentialContourPlot(X,Y,Vfree,10) %#ok<UNRCH>
end

%% Initial State

% ground state of SHO using Pethick & Smith, eq. 6.17
% (for GPE, use Thomas-Fermi distribution instead)
ax          = sqrt(hbar/(m*omegax))*kL;     % oscillator lengths
ay          = sqrt(hbar/(m*omegay))*kL;
amp         = 1/(pi*ax*ay)^(1/2);           % amplitude for normalization

% initial wave function matrix (normalized)
cii         = amp.*exp(-1/2.*((X./ax).^2.+(Y./ay).^2));

c           = zeros(Nx,Ny,Rsaves+Fsaves);	% matrix to hold saved wave functions
c(:,:,1)    = cii;                          % store initial wave function
% note Rsave == 1, so this first save is already accounted for

% plot of initial state
if makeInitialStatePlot
    makePositionPlot(X,Y,cii,xMax,yMax) %#ok<UNRCH>
end

%% Propagators

% for interactions, redefine UV at every step
UVramp  = @(t) exp(-1i*Vramp(t)*dtR/2);	% half position propagator (B/2) for ramp
UVfree  = exp(-1i*Vfree*dtR/2);         % half position propagator (B/2) for free evolution
k       = sqrt(Kx.^2+Ky.^2);
Uk      = exp(-1i*k.^2*dtF);            % momentum propagator (A)
UkShf   = fftshift(Uk);                 % shifted momentum propagator to account for shift in momentum matrix

%% Time Evolution

tic
t = 0;

% ramp
for ii = 1:tRsteps
    cii = UVramp(t).*cii;           % propagate in position
    cii = ifft2(UkShf.*fft2(cii));	% transform to momentum space, propagate, transform back
    cii = UVramp(t).*cii;           % propagate in position
    if ii/Rsavestep >= Rsave
        fprintf('evolving ramp step %d of %d\n',ii,tRsteps);
        Rsave = Rsave+1;
        c(:,:,Rsave) = cii;          % save wave function
    end
    t = t+dtR;
end

disp(Rsave);

% free evolution
for ii = 1:tFsteps
    cii = UVfree.*cii;              % propagate in position
    cii = ifft2(UkShf.*fft2(cii));  % transform to momentum space, propagate, transform back
    cii = UVfree.*cii;              % propagate in position
    if ii/Fsavestep >= Fsave
        fprintf('evolving free evolution step %d of %d\n',ii,tFsteps);
        Fsave = Fsave+1;
        if Fsave <= Fsaves          % this is hard-coded to fix an error where Fsave goes over by 1 sometimes
            c(:,:,Rsaves+Fsave) = cii; % save wave function
        end
    end
    t = t+dtF;                      % not actually used here, but will need for interactions
end
toc

%% Make a gif
tic

% initialize relevant objects
fGIF    = figure;
axGIF   = axes;
im      = cell(1,Rsaves+Fsaves);

% output file name in form yymmdd-HHMM_Nx[]_dx[]_dt[]_tf[].gif
filename = sprintf('%s_%s%g_%s%g_%s%g_%s%g.gif',datestr(now,'yymmdd-HHMM'),'Nx',Nx,'dx',dx,'dt',dtF,'tf',tFree);

speedGIF = 0.01;                          % playback speed of GIF; 1 == true speed
delayR   = tRamp_ms/1000/Rsaves/speedGIF;	% make the gif at [speed] x the actual rate
delayF   = tFree_ms/1000/Fsaves/speedGIF;

% annotate plot with parameters
params1	= {['Nx: ' num2str(Nx)];
           ['dx: ' num2str(dx)];
           ['ax: ' num2str(ax)];
           ['ay: ' num2str(ay)]};     
params2 = {['dtR:   ' num2str(dtR_ms) ' ms'  ];
           ['dtF:   ' num2str(dtF_ms) ' ms'  ];
           ['tRamp: ' num2str(tRamp_ms) ' ms'];
           ['tFree: ' num2str(tFree_ms) ' ms']};
params3 = {['Fx: ' num2str(FX)             ];
           ['Fy: ' num2str(FY)             ];
           ['speed: ' num2str(speedGIF) 'x']};
annotation('textbox',[.10 .88 .2 .1],'String',params1,'EdgeColor','none','FontName','FixedWidth');
annotation('textbox',[.25 .88 .2 .1],'String',params2,'EdgeColor','none','FontName','FixedWidth');
annotation('textbox',[.45 .88 .2 .1],'String',params3,'EdgeColor','none','FontName','FixedWidth');

% plot all saved wavefunctions and render to image
for ii = 1:(Rsaves+Fsaves)
    fprintf('rendering frame %d\n',ii);
    timeAnnotation(ii,Rsaves,Fsaves,tRamp_ms,tFree_ms,0.1);
    cla;                        % otherwise plots are stacked on top of each other
    set(fGIF,'color','w');      % otherwise you get a gray background in the pdf/gif
    hold on
    sSaves = surf(X,Y,abs(c(:,:,ii)).^2);
    xlim([-xMax xMax]);
    ylim([-yMax yMax]);
    zlim([0 3*amp^2]);
    labelparams = {'interpreter','latex','FontSize',24};
    xlabel('$x\ (1/k_L)$','Rotation',17,labelparams{:});
    ylabel('$y\ (1/k_L)$','Rotation',-56,labelparams{:});
    zlabel('$c^2$',labelparams{:})
    sSaves.EdgeColor = 'none';
    view([-25.5 75.1]);
    contour(X,Y,Vlat,[0 5 10 14],"LineColor",[1,0,0]);
    axGIF.CLim = [0,.5*amp^2];     % define color scale
    hold off
    frame   = getframe(fGIF);
    im{ii}  = frame2im(frame);
end
toc

% append images to GIF file
tic

for ii = 1:(Rsaves+Fsaves)
    fprintf('appending image %d\n',ii);
    [A,map] = rgb2ind(im{ii},256);
    if ii == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delayR);
    else
        if ii <= Rsaves
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delayR);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delayF);
        end
    end
end
toc

%% Functions

function makePotentialSurfacePlot(X,Y,V,xMax,yMax) %#ok<DEFNU>
    f1 = figure;
    s1 = surf(X,Y,V);
    xlim([-xMax xMax]);
    ylim([-yMax yMax]);
    labelparams = {'interpreter','latex','FontSize',24};
    xlabel('$x\ (1/k_L)$','Rotation',17,labelparams{:});
    ylabel('$y\ (1/k_L)$','Rotation',-56,labelparams{:});
    zlabel('$V(x,y)\ (E_r)$',labelparams{:})
    s1.EdgeColor = 'none';
    view([-25.5 75.1])
end

function makePotentialContourPlot(X,Y,V,Nlevels) %#ok<DEFNU>
    f2 = figure;
    contourf(X,Y,V,Nlevels)
    colorbar
    colormap(parula)
    set(f2,'Units','Pixels');
    labelparams = {'interpreter','latex','FontSize',24};
    xlabel('$x\ (1/k_L)$',labelparams{:});
    ylabel('$y\ (1/k_L)$',labelparams{:});
    ylabel(colorbar,'$V(x,y)\ (E_r)$',labelparams{:});
end

function makePositionPlot(X,Y,cii,xMax,yMax) %#ok<DEFNU>
    f3 = figure;
    ax3 = axes;
    s3 = surf(X,Y,cii);
    xlim([-xMax xMax]);
    ylim([-yMax yMax]);
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    zlabel('$c$','interpreter','latex')
    s3.EdgeColor = 'none';
    view([-25.5 75.1])
end

function timeAnnotation(ii,Rsaves,Fsaves,tRamp_ms,tFree_ms,unitTime_ms)
    % t_final should be in ms
    % unitTime_ms is the interval in ms at which you want the diplayed time
    % on the GIF to increase; only multiples of this value will be displayed
    if ii <= Rsaves
        t_per_frame = tRamp_ms/Rsaves;
        displayTime = round(ii*t_per_frame/unitTime_ms)*unitTime_ms;
        annotation('textbox',[.2 .7 .15 .08],'String',{[num2str(displayTime) ' ms'];'ramp'},'EdgeColor','none','FontName','FixedWidth','BackgroundColor','white');
    else
        t_per_frame = tFree_ms/Fsaves;
        displayTime = round((ii-Rsaves)*t_per_frame/unitTime_ms)*unitTime_ms + tRamp_ms;
        annotation('textbox',[.2 .7 .15 .08],'String',{[num2str(displayTime) ' ms'];'free evolution'},'EdgeColor','none','FontName','FixedWidth','BackgroundColor','white');
    end
end