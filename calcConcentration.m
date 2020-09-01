function [C, rs, ETA, xcIter, cIter, xcFinal, beta, skin, kDaverage, kDamage, t ] = ...
    calcConcentration ( k, dP, phi, kc, L, Sor, phic, rhoc, Csolid, muF, Cf,...
    kTao, kPrime, FI, taoCrit, f, n, dx, dt, time, x, CdStar,...
    A, measuredFL, tMax, cn, betaInput, gSwitch)
cIter = 1;
kinMilidarcy = k /(9.869233e-16);
kcinMilidarcy = kc / (9.869233e-16);
Swi = 100*phi^2/(100*phi^2 + sqrt(kinMilidarcy));  % Coates and Denoo's permeability model
waterFraction = phi*(1 - Swi - Sor);
rhoF = 1000; % kg/m3, assumed water density

% DETERMINE INERTIAL COEFFICIENT
if strcmp(betaInput, 'liu') == 1
    beta = 2.881816036000000e-06/k/phi;   % 1/m,  Liu et al. (1995) ref: Civan
    betac = 2.881816036000000e-06/kc/phic; % 1/m
elseif strcmp(betaInput, 'FV') == 1
    beta = 4.1e11/kinMilidarcy^(1.5);     % 1/m, Friedel and Voigt (2006)
    betac = 4.1e11/kcinMilidarcy^(1.5);
else
    beta = betaInput;
    betac = betaInput;
end
% termND = rhoF*beta*dP/L;

% find g index for the diffusion term
if strcmp(gSwitch, 'off') == 1
    g = findG ( kc * 1e12 );
else
    g = gSwitch;
end

% Concentration I.C. and B.C.s
C = zeros(n, 1);
C(1) = Cf;
FL = 0;

% Time loop
t = 1;
xcTol = 1e-11;
uOld = 0;
xc = 0;
while FL <= measuredFL
    % Velocity loop
    xcError = inf;
    xcIter = 0;
    xc(t+1) = 0;
    while xcError > xcTol  % this loop outputs the velocity
        uRoots = getVelocity ( kc, k, xc(t+1), L, dP, muF, rhoF, phi, phic, beta, betac  );
        if length(uRoots) == 2
            u = uRoots(2);
        else
            u = uRoots;
        end
        tao = findTao ( u*100, kPrime, FI, taoCrit);
%         tao = findTao ( 125, kPrime, FI, taoCrit);
        xcOld = xc(t+1);
        xc(t+1) = xc(t) + (u*(1-phi)*Csolid-kTao*tao*1e-1) * dt / ( (1-phic)*rhoc );
%         xc(t+1) = xc(t) + (u*Csolid-kTao*tao*1e-1) * dt / ( (1-phic)*rhoc );
        xcError = abs(xc(t+1) - xcOld);
        xcIter = xcIter + 1;
        if xcIter > 50000
            error('Velocity and Thickness could not be coupled!');
            return
        end
    end
    
    % FX of pressure on cake compaction
%     [phic, kc] = compressCake (u, muF, xc(t+1), kc, phic ); 
    
    % calculate the current fluid loss volume ( trapezoidal integration )
    FL = dt * (u + uOld)/2 * A + FL;
    uOld = u;
    
    % Skin Calculation
    [skin(t), kDaverage, kDamage] = findSkin (C, CdStar, Cf, Sor, Swi, x, n, k );
    
    % % Concentration Calculation % %
    diffusionCoeff = (f*u^g)*dt/dx^2; % diffusion term coefficient
    convectionCoeff = u*dt/2/dx/waterFraction; % convection term coefficient
    % coefficients matrix diagonals
    up = (convectionCoeff-diffusionCoeff);
    mid = (1+2*diffusionCoeff);
    down = -(diffusionCoeff+convectionCoeff);
    % solving the diffusion convection equation
    bVector = C(2:n-1);
    bVector(1) = C(2) + (diffusionCoeff+convectionCoeff)*C(1);
    bVector(n-2) = C(n-1) - (convectionCoeff-diffusionCoeff)*C(n-2);
    C(2:n-1) = Thomas(down*ones(1, n-3), mid*ones(1, n-2), up*ones(1, n-3), bVector);
%         [C(2:n-1), cIter] = calcTriGaussSeidel (down, mid, up, bVector, zeros(n-2, 1), 1e-6, 1);    

    % Flow into the core animation
%         figure(1)
%         set(gcf,'color','w');
%         set(gcf,'Position',[100 100 475 700]);
%         plot(x./L, C./Cf, 'r' ,'LineWidth', 2);
%         title('Dimensionless Mud Filtrate Concentration');
%         grid on
%         xlabel('Dimensionless Position');
%         ylabel('Dimensionless Concentration');
%         ylim([0 1]);
%         text(.62, .95, ['Elapsed time: ', num2str(time(t)), ' s']);

    if time(t) >= tMax
        error('Maximum allowable time reached!');
        return
    end
    t = t + 1;
end


% finding invasion depth = max{x: C > Cd}
rs = 0;
for i = 1:n
    if C(i) > CdStar*Cf
        rsOld = rs;
        rs = x(i);
        if rs >= rsOld
            rsFinal = rs;
        end
    end
    
end
rs = rsFinal*100; % cm
ETA = time(t);
xcFinal = xc(t);
end

function taoPrime = findTao ( v, kPrime, FI, taoCrit)
% Reservoir Formation Damage, F. Civan (2015) p. 299 & p.309 Table 12.1
tao = kPrime*(8*v)^FI;
if tao < taoCrit
    taoPrime = 0;
else
    taoPrime = tao - taoCrit;
end
end

function g = findG ( kc )
% kc is in micron2

if kc <= 2e-9
    g = 1.47;
elseif kc > 2e-9 && kc <= 1e-8
    g = 1.47 + .0225*(kc*1e9 - 2);
elseif kc > 1e-8 && kc <= 2e-8
    g = 1.65 + .03*(kc*1e9 - 10);
elseif kc > 2e-8
    g = 1.95;
    
end
end

%% PROJECT EXTENSIONS

function uRoots = getVelocity (kc, k, xc, L, dP, muf, rhoF, phi, phic, beta, betac )

%         termD = 1/ ( k*dP/muf/L * (kc/k)/(kc/k + xc/L) );
%         uNondarcy = [termND termD -1];
        averagePorosity = 1 - (((1-phi)*L + (1-phic)*xc)/ (xc + L));
        tortuosity = 1.16*averagePorosity^1.15 ; % (Winsauer et al., 1952)
        termU = muf*(xc/kc + L/k);
%         termU2 = 2.881816036e-6*rhoF*tortuosity/averagePorosity*(xc/kc+L/k);
        termU2 = (1+betac*xc/L/(beta))*rhoF*beta*tortuosity*L;
        uNondarcy = [termU2 termU -dP];
        uRoots = roots(uNondarcy);
end

function [skin, kDamageAverage, kDamage] = findSkin (C, CdStar, Cf, Sor, Swi, x, n, k )

rs = 0;
for i = 1:n
    if C(i) > CdStar*Cf
        rsOld = rs;
        rs = x(i);
        if rs >= rsOld
            rsFinal = rs;
        end
    end
end

no = 2;
for i = 2:n
    % Oil Relative Permeability calculation
    SwMax = 1 - Sor;
    SwMin = Swi;
    Smud(i) = C(i)/Cf*(SwMax - SwMin) + SwMin;
%     kro(i) = kroMax *   ( (1 - Smud(i) - Sor) / ( 1 - Swi - Sor ) )^no;
    kDamage(i) = k *   ( (1 - Smud(i) - Sor) / ( 1 - Swi - Sor ) )^no;
end

% Skin Calculation
permSum = 0;
for j = 2:n
    if x(j) <= rsFinal
        permSum = permSum + (kDamage(j)+kDamage(j-1))/2*(x(j)-x(j-1));
    end
end
kDamageAverage = ( permSum ) / ( rsFinal);
skin = k/kDamageAverage*rsFinal - rsFinal;
end

function [newPhic, newKc, Pc] = compressCake (u, muF, xc, kc, phic )

Pc = u*muF*xc/kc; % differential pressure inside the cake
% delta = 0.1;
% nu = 0.4;
% newPhic = phic / Pc^(delta*nu);
% newKc = kc / Pc^nu;
c = 5e-6 * .00064516/4.44822; % formation compressiblity, m2/N
newPhic = phic - c*phic*Pc;
r = .6 * 1e-9; % nanometers to m (0.6 for first set, 2 for second set)
% r = 0.0314*(kc /(9.869233e-16)/newPhic)^.5;
Sp = 2/r;
Sg = Sp*(newPhic/(1-newPhic));
newKc = 1/5/Sg^2 * (newPhic^3/ (1-newPhic)^2); % Carman-Kozeny relationship

end