function PISA_calibration_plots(WTG)
%close all;
for loc = 1:length(WTG)
 
%% Triple plot 
excelname1 = '\\cowi.net\projects\A125000\A127858\GEO\PISA_SoilSprings\CalibrationOfPISASoilReactionCurves\Validation_Homogeneous_soil\COSPIN_results.xlsx';
[NUM] = xlsread(excelname1,WTG{loc});
excelname2 = '\\cowi.net\projects\A125000\A127858\GEO\PISA_SoilSprings\CalibrationOfPISASoilReactionCurves\Validation_Homogeneous_soil\StructuralForces_Plaxis.xlsx';
[PLAX] = xlsread(excelname2,['' WTG{loc} '_D200' '']);
excelname3 = '\\cowi.net\projects\A125000\A127858\GEO\PISA_SoilSprings\CalibrationOfPISASoilReactionCurves\Validation_Homogeneous_soil\FxUx_3DFE_D10.xlsx';
[PLAX_fxux] = xlsread(excelname3,['' WTG{loc} '_mudl' '']);

figure('visible','on', 'Position', [10 10 900 600])
subplot(1,3,1)
hold on
plot(NUM(:,1),NUM(:,4))
plot(PLAX(:,4)*2, PLAX(:,1), '-rx')
ylabel('Depth [m]')
xlabel('Moment [kNm]')
grid on
legend('COSPIN - 1D', 'PLAXIS - 3D','Location','southeast')

subplot(1,3,2)
hold on
plot(NUM(:,2),NUM(:,4))
plot(PLAX(:,3)*2, PLAX(:,1), '-rx')
ylabel('Depth [m]')
xlabel('Shear [kNm]')
grid on
legend('COSPIN - 1D', 'PLAXIS - 3D','Location','northeast')

subplot(1,3,3)
hold on
plot(NUM(:,3),NUM(:,4))
plot(PLAX(:,11)*0.5, PLAX(:,10)*-1, '-rx')
ylabel('Depth [m]')
xlabel('Displacement [mm]')
grid on
legend('COSPIN - 1D', 'PLAXIS - 3D','Location','southeast')

savefig(['' '\\cowi.net\projects\A125000\A127858\GEO\Spring_Calibration\Upper_clays\results\moment_shear_disp_plots\' WTG{loc} 'mom_shea_disp_plot' ''])
saveas(gcf, ['' '\\cowi.net\projects\A125000\A127858\GEO\Spring_Calibration\Upper_clays\results\moment_shear_disp_plots\' WTG{loc} 'mom_shea_disp_plot.png' ''])

%% Load-Displacement plot 1 - limit 2000mm

figure('visible','on')
hold on
plot(NUM(:,6),NUM(:,7))
plot(PLAX_fxux(:,5)*1000, PLAX_fxux(:,6)*2, '-rx')
grid on
eval(['title(''Pile head load-displacement curve plot, ' WTG{loc} ' '');'])
xlim([0 2000])
legend('COSPIN - 1D', 'PLAXIS - 3D','Location','southeast')
ylabel('Pile head load ratio H [-]')
xlabel('Pile head deflection [mm]')

savefig(['' '\\cowi.net\projects\A125000\A127858\GEO\Spring_Calibration\Upper_clays\results\load_disp_plots\0-2000mm\' WTG{loc} 'load_disp_plot_2000' ''])
saveas(gcf, ['' '\\cowi.net\projects\A125000\A127858\GEO\Spring_Calibration\Upper_clays\results\load_disp_plots\0-2000mm\' WTG{loc} 'mom_shea_disp_plot_2000.png' ''])

%% Load-Displacement plot 2 - limit 20mm

figure('visible','on')
hold on
plot(NUM(:,6),NUM(:,7))
plot(PLAX_fxux(:,5)*1000, PLAX_fxux(:,6)*2, '-rx')
grid on
eval(['title(''Pile head load-displacement curve plot, ' WTG{loc} ' '');'])
xlim([0 20])
legend('COSPIN - 1D', 'PLAXIS - 3D','Location','southeast')
ylabel('Pile head load ratio H [-]')
xlabel('Pile head deflection [mm]')

savefig(['' '\\cowi.net\projects\A125000\A127858\GEO\Spring_Calibration\Upper_clays\results\load_disp_plots\0-20mm\' WTG{loc} 'load_disp_plot_20' ''])
saveas(gcf, ['' '\\cowi.net\projects\A125000\A127858\GEO\Spring_Calibration\Upper_clays\results\load_disp_plots\0-20mm\' WTG{loc} 'mom_shea_disp_plot_20.png' ''])

end

end