figure(100)
plot(Es{1,1}(:,2,2),node.level(2:end),'k',Es{1,1}(:,1,2),node.level(1:end-1),'r')
title('Shear force')
xlabel('Shear force [kN]')
ylabel('Depth [m]')
legend('bottom of element','top of element')
grid on

figure(101)
plot(Es{1,1}(:,2,3),node.level(2:end),'k',Es{1,1}(:,1,3),node.level(1:end-1),'r')
title('Bending moment')
xlabel('Bending moment [kNm]')
ylabel('Depth [m]')
legend('bottom of element','top of element')
grid on

figure(102)
plot(output.m_UR(:,2),node.level(2:end-1),'k',output.m_UR(:,1),node.level(1:end-2),'k')
title('distributed moment')
xlabel('Distributed moment [kNm/m]')
ylabel('Depth [m]')
legend('bottom of element','top of element')
grid on

figure(103)
plot(Es{1,1}(1:end-1,2,2)-(output.m_UR(:,2)+output.m_UR(:,1))/2,node.level(2:end-1),'k',Es{1,1}(1:end-1,1,2)-(output.m_UR(:,2)+output.m_UR(:,1))/2,node.level(1:end-2),'r')
title('Shear force, corrected')
xlabel('Shear force [kN]')
ylabel('Depth [m]')
legend('bottom of element','top of element')
grid on

figure(104)
plot(Es{1,1}(1:end-1,2,2)-(output.m_UR(:,2)+output.m_UR(:,1))/2,node.level(2:end-1),'k',Es{1,1}(1:end-1,1,2)-(output.m_UR(:,2)+output.m_UR(:,1))/2,node.level(1:end-2),'r',Es{1,1}(:,2,2),node.level(2:end),'k:',Es{1,1}(:,1,2),node.level(1:end-1),'r:')
title('Shear force, corrected vs uncorrected')
xlabel('Shear force [kN]')
ylabel('Depth [m]')
legend('bottom of element - corrected','top of element - corrected','bottom of element - uncorrected','top of element - uncorrected')
grid on