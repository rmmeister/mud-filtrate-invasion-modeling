%% INITIALIZE MATLAB
clear all
clc
close all
format long

%% DASHBOARD
% Initial rock and fluid data are entered in this section.

% Table 1
overbalancePressure = [2.76 1.50 2.76 2.76 2.41 3.90 2.10 1.50 2.76 1.50...
    3.10 3.10 3.10 2.76] .* 1e6;                        % Pa
measuredFluidLoss = [2.33 7.37 1.77 2.53 5.17 2.44 7.98 5.34 2.78 7.51 ...
    2.52 2.29 2.61 2.78] .* 1e-6;                       % m3
corePermeability = [33.09 64.81 18.96 170.69 128.23 153.27 72.03 45.32 ...
    329.77 113.29 710.21 379.43 307.45 595.61] * 1e-15; % 1e-3 micron2 to m2
corePorosity = [.075 0.0953 0.1451 0.2482 0.2376 0.2434 0.2240 0.2440 ...
    0.2953 0.0932 0.3025 0.2804 0.2895 0.2965];         % fraction
cakePermeability = [2.220e-9 3.960e-8 1.180e-9 2.120e-9 1.026e-8 ...
    1.400e-9 2.852e-8 1.741e-8 2.390e-9 4.135e-8 1.741e-9 1.477e-9 ...
    1.895e-9 2.390e-9 ] * 1e-12;                        % micron2 to m2
measuredInvasionDepth = [6.62 6.50 5.72 5.66 5.53 6.31 5.87 4.85 5.62 ...
    6.61 5.45 5.54 5.91 5.57];                          % cm
calculatedInvasionDepth = [6.95 6.49 5.39 5.17 5.91 6.16 6.26 4.61 5.09 ...
    6.57 5.32 5.22 5.345 5.148];                        % cm

% TABLE 2
Length = 4 / 39.37007874;        % in to meters
Area = 20.42 * 1e-4;             % cm2 to m2
residualOilSaturation = .3;      % fraction
cakePorosity = .1;               % fraction
cakeDensity = 2440;              % kg/m3
solidConcentration = 40;         % kg/m3
filtrateViscosity = 1 * 1e-3;    % kg/ms
filtrateConcentration = 1000;    % kg/m3
CdStar = .001;                   % dimensionless

% OTHER
kTao = 0 * 100;             % s/cm to s/m
kPrime = 8;               % dyne/cm2/s^FI to N/m2/s^FI * 1e-1
FI = .319;                       % dimensionless
criticalShearStress = 5;  % dyne/cm2 to N/m2 * 1e-1
f = 51.7;
intertialCoeff = 'FV'; % 1/m, put down 'liu' to calculate using Liu et al. (1995) formula
gSwitch = 'off';
%% OPERATION SELECTION
% This section is about picking the desired operation.
while true
close all
clc
% Choose Operation
Operation = input(['Enter:\n(1) for Single-Core Simulation,\n(2) ', ...
    'for Error Analysis,\n(3) for Sensitivity Analysis,\nPress (enter) to quit.\n-> ']);
if Operation == 1
    coreOption = input('Enter Core no. (1 to 14): ');
elseif Operation == 2
    coreOption = 1:1:14;
elseif Operation == 3
    coreOption = input('Enter Core no. (1 to 14): ');
    sensOption = input(['\nSelect a parameter to perform sensitivity test on:\n', ...
          '(a):    Area           |', '(k):    Core Permeability' ...
        '\n(v):    Viscosity      |', '(p):    Overbalance Pressure',...
        '\n(phi):  Core Porosity  |', '(cs):   Solid Concentration',...
        '\n(phic): Cake Porosity  |', '(rhoc): Cake Density', ...
        '\n(f):    f (of D)       |', '(g):    g (of D)', '\n-> '], 's');
    if strcmp(sensOption, 'a') == 1
        sensParameter = input('Insert "Cross-Sectional Area" values to compare (cm2): ');
    elseif strcmp(sensOption, 'v') == 1
        sensParameter = input('Insert "Viscosity" values to compare (cp): ');
    elseif strcmp(sensOption, 'phi') == 1
        sensParameter = input('Insert "Core Porosity" values to compare (fraction): ');
    elseif strcmp(sensOption, 'phic') == 1
        sensParameter = input('Insert "Cake Porosity" values to compare (fraction): ');
    elseif strcmp(sensOption, 'k') == 1
        sensParameter = input('Insert "Core Permeability" values to compare (mD): ');
    elseif strcmp(sensOption, 'p') == 1
        sensParameter = input('Insert "Over-balance Pressure" values to compare (MPa): ');
    elseif strcmp(sensOption, 'rhoc') == 1
        sensParameter = input('Insert "Cake Density" values to compare (kg/m3): ');
    elseif strcmp(sensOption, 'cs') == 1
        sensParameter = input('Insert "Solid Concentration" values to compare (kg/m3): ');
    elseif strcmp(sensOption, 'f') == 1
        sensParameter = input('Insert "f" values to compare (m): ');
    elseif strcmp(sensOption, 'g') == 1
        sensParameter = input('Insert "g" values to compare: ');
    else
        fprintf('\nEnter only letters written inside the parantheses!\n');
    end
    lenS = length(sensParameter);
else
    break
end

%% SIMULATION
% Simulation of mud invasion into the desired core/cores.

% Meshing Option
gridNo = 1000; % number of nodes
dx = Length/gridNo;
x = dx/2:dx:Length-dx/2;
dt = 5; % sec, 0.05 is OK for min errors
tMax = 40000; % sec
t = 0:dt:tMax;
xStar = x./Length;

% Run simulation
if Operation ~= 3
    fprintf('Running Simulation... \n\n' );
    for coreNo = coreOption
        [C, predictedInvasionDepth(coreNo), ETA, xcIter, cIter, xc, beta, skin, averageDamagePerm, kDamage,...
            tEnd ] =  calcConcentration...
            ( corePermeability(coreNo), overbalancePressure(coreNo), corePorosity(coreNo), ...
            cakePermeability(coreNo), Length, residualOilSaturation, cakePorosity, ...
            cakeDensity, solidConcentration, filtrateViscosity, ...
            filtrateConcentration, kTao, kPrime, FI, criticalShearStress, ...
            f, gridNo, dx, dt, t, x, CdStar, Area, measuredFluidLoss(coreNo), tMax, coreNo, ...
            intertialCoeff, gSwitch);

        invasionError(coreNo, 1) = abs(predictedInvasionDepth(coreNo) - ...
            measuredInvasionDepth(coreNo))/measuredInvasionDepth(coreNo)*100;
    end
    fprintf('Simulation ended successfully.\n\n');
end
%% SINGLE-CORE SIMULATION RESULTS
% Error is obtained and the concentration profile is generated for chosen
% core sample.
if Operation == 1
    fprintf('The results for Core no.(%d) are as follows:\n', coreNo);
    fprintf('Predicted Invasion Depth: %g cm \n', predictedInvasionDepth(coreNo));
    fprintf('Measured Invasion Depth:  %g cm \n', measuredInvasionDepth(coreNo));
    fprintf('Invasion Depth Error:     %g percent \n', invasionError(coreNo, 1));
    fprintf('Injection Time:           %g minutes \n', ETA/60 );
    fprintf('Skin:                     %g cm \n', skin(tEnd-1)*100);
    fprintf('Damage zone permeability: %g mD\n\n', averageDamagePerm /(9.869233e-16));
    % C vs. x plot setting
    pltOption = input(['Concentration Profile:       (1)\n',...
                       'Damage Permeability Profile: (2)\n', ...
                       'Skin vs. Elapsed Time:       (3)\n', '->']);
    if pltOption == 1
        figure(1)
        set(gcf,'Position',[300 300 475 700]);
        plot(xStar, C./filtrateConcentration, 'k' ,'LineWidth', 2);
        grid on
        xlabel('Dimensionless Position');
        ylabel('Dimensionless Concentration');
        title('Dimensionless Mud Filrate Concentration');
        ylim([0 1]);
        legend(['Core no. ', num2str(coreNo)]);
        fprintf('Press any key to move on');
        pause
    elseif pltOption == 2
        figure
        plot(x.*100, kDamage /(9.869233e-16), 'LineWidth', 2);
        xlabel('Core Length, cm');
        ylabel('Damage Zone Permeability, mD');
        xlim([0, predictedInvasionDepth(coreNo)]);
        grid on;
        fprintf('Press any key to move on');
        pause
    elseif pltOption == 3
        figure
        plot(t(1:tEnd-1)./60, skin*100, 'r', 'LineWidth', 2);
        xlabel('Elapsed Time, min');
        ylabel('Skin, cm');
        grid on;
        fprintf('Press any key to move on');
        pause
    end
end
%% INVASION DEPTH COMPARISON
% Error Analysis of all core samples and comparison to experimental data.
if Operation == 2
    % invasion error
    errPlot = figure('Position', [400 400 700 200]);
    axes1 = axes('Parent',errPlot,'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14]);
    box(axes1,'on');
    hold(axes1,'all');
    plot(1:1:14, invasionError, 'LineWidth', 2, 'Marker', '^', 'Color', 'r'...
        , 'LineStyle', 'none');
    grid on;
    title('Predicted Invasion Depth Error Relative to Measured Invasion Depth');
    xlabel('Core No.');
    ylabel('Invasion Depth Relative Error, percent');

    % invasion depth comparison
    figure
    plot(linspace(4.5, 7.1), linspace(4.5, 7.1), 'k', 'LineWidth', 2);
    hold on
    plot(measuredInvasionDepth, calculatedInvasionDepth, 'k', 'LineWidth',...
        1,'Marker', 'square', 'MarkerFaceColor',[0 0 0], 'LineStyle', 'none');
    plot(measuredInvasionDepth, predictedInvasionDepth, 'k^', 'LineWidth', 1,...
        'MarkerFaceColor',[0 0 0], 'LineStyle', 'none');
    hold off
    xlim([4.5 7.1]);
    ylim([4.5 7.1]);
    legend('Measured', 'Empirical', 'Model', 'Location', 'NorthWest');
    xlabel('Measured');
    ylabel('Predicted');

    % invasion table
    invasionTableOpt = input(['Display the Model Invasion Depth table?'...
        ,'\n(y): yes\n(s): save data as Excel spreadsheet\n(enter): no and move on.\n-> '], 's');
    VarNames = {'CoreNo', 'PredictedInvasionDepth', 'RelativeError', ...
        'CoreNo_contd', 'PredictedInvasionDepth_contd', 'RelativeError_contd' };
    VarUnits = {'', 'cm', '%', '', 'cm', '%'};
    format bank
    invasionTable = table(coreOption(1:7)', predictedInvasionDepth(1:7)',...
        invasionError(1:7), coreOption(8:14)', predictedInvasionDepth(8:14)'...
        , invasionError(8:14));
    invasionTable.Properties.VariableNames = VarNames;
    invasionTable.Properties.VariableUnits = VarUnits;
    if strcmp(invasionTableOpt, 'y') == 1 || strcmp(invasionTableOpt, 'Y') == 1
        disp(invasionTable);
        fprintf('Press (Enter) to move on');
        pause
    elseif strcmp(invasionTableOpt, 's') == 1 || strcmp(invasionTableOpt, 'S') == 1
        tableXL = 'invasionTable.xlsx';
        writetable(invasionTable, tableXL, 'Sheet', 1);
    end
    format long
end
%% SENSITIVITY TEST
% Sensitivity Analysis on chosen parameters on the desired core sample.
if Operation == 3
    fprintf('\nRunning Sensitivity... \n\n' );
    for coreNo = coreOption
        sensPlot = zeros(length(lenS));
        for i = 1:lenS
            if strcmp(sensOption, 'a') == 1
                Area = sensParameter(i) * 1e-4;
                sensName = 'Cross-sectional Area';
            elseif strcmp(sensOption, 'v') == 1
                filtrateViscosity = sensParameter(i) * 1e-3;
                sensName = 'Filtrate Viscosity';
            elseif strcmp(sensOption, 'phi') == 1
                corePorosity(coreNo) = sensParameter(i);
                sensName = 'Core Porosity';
            elseif strcmp(sensOption, 'phic') == 1
                cakePorosity = sensParameter(i);
                sensName = 'Cake Porosity';
            elseif strcmp(sensOption, 'k') == 1
                corePermeability(coreNo) = sensParameter(i) * 1e-15;
                sensName = 'Core Permeability';
            elseif strcmp(sensOption, 'p') == 1
                overbalancePressure(coreNo) = sensParameter(i) * 1e6;
                sensName = 'Overbalance Pressure';            
            elseif strcmp(sensOption, 'rhoc') == 1
                cakeDensity = sensParameter(i);
                sensName = 'Cake Density';
            elseif strcmp(sensOption, 'cs') == 1
                solidConcentration = sensParameter(i);
                sensName = 'Solid Concentration';
            elseif strcmp(sensOption, 'f') == 1
                f = sensParameter(i);
                sensName = 'f (of D=fu^g)';
            elseif strcmp(sensOption, 'g') == 1
                gSwitch = sensParameter(i);
                sensName = 'g (of D=fu^g)';
            end
            % run for selected parameter value
            [C, predictedInvasionDepth(coreNo), ETA(i), xcIter, cIter, xc(i) ] =  calcConcentration...
                ( corePermeability(coreNo), overbalancePressure(coreNo), corePorosity(coreNo), ...
                cakePermeability(coreNo), Length, residualOilSaturation, cakePorosity, ...
                cakeDensity, solidConcentration, filtrateViscosity, ...
                filtrateConcentration, kTao, kPrime, FI, criticalShearStress, ...
                f, gridNo, dx, dt, t, x, CdStar, Area, measuredFluidLoss(coreNo), tMax, coreNo,...
                intertialCoeff, gSwitch);
            figure(1)
            sensColor = {'k','b','r','g','y','c','m',[.5 .6 .7],[.8 .2 .6], [.8 .2 .5]};
            set(gcf,'Position',[300 300 700 700]);
            sensPlot(i) = plot(xStar, C./filtrateConcentration, 'k' ,'LineWidth', 2);
            grid on
            xlabel('Dimensionless Position');
            ylabel('Dimensionless Concentration');
            title('Dimensionless Mud Filrate Concentration');
            ylim([0 1]);
            hold on
            set(sensPlot(i), 'DisplayName', [sensName, ' = ', num2str(sensParameter(i)),...
                ', t_{inj} = ', num2str(round(ETA(i)/60)), ' min'], 'Color', sensColor{i});
        end
    end
    legend('show');
    fprintf('Press any key to exit\n');
    pause
    break
end
end