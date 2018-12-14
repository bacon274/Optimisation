%This script is used to find MPP

pmax = max(PV.signals.values(:,2));
vrange = max(PV.signals.values(:,1));
irange = max(PV.signals.values(:,3));
[tf,index]=ismember(pmax,PV.signals.values(:,2));
disp(' MPP power: ')
disp(PV.signals.values(index,2));
disp(' MPP voltage: ')
disp(PV.signals.values(index,1));
disp(' MPP current: ');
disp(PV.signals.values(index,3));
figure(1)
plot(PV.signals.values(:,1),PV.signals.values(:,2)); 
axis([0 vrange 0 pmax]);
figure(2)
plot(PV.signals.values(:,1),PV.signals.values(:,3));
axis([0 vrange 0 irange]);