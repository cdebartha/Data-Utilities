a = load('bp1p8_flowRate.dat');
m_x1 = a(:,1);
m_x2 = a(:,2);
m_x3 = a(:,3);
t = linspace(0.0,11.0,12);

%%
figure(1)
plot(t, m_x1, '-*', 'Linewidth', 1.2);
hold on
plot(t, m_x2, '-^', 'Linewidth', 1.2);
hold on
plot(t, m_x3, '-s', 'Linewidth', 1.2);
grid on
legend('$x = 0.0572$, intake inlet', '$x = 1.0$, merger section', '$x = 2.8$, intake outlet','Interpreter', 'latex',...
    'location','southeast');
legend boxoff;
ylim([0.056, 0.065]);
xlim([0, 12])
% xlim([5.0, 11.0]);

ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
ax.TickLabelInterpreter = 'latex';
set(gcf, 'Position',  [100, 100, 875, 350]);
%title('Mass flow rate variation with time at three different locations for $p_b^* = 1.4$','Interpreter', 'latex');
xlabel('t', 'Interpreter', 'latex');
ylabel('$\dot{m}$', 'Interpreter', 'latex');
%%
a1 = load('bp1p4_flowRate.dat');
a2 = load('bp1p5_flowRate.dat');
a3 = load('bp1p6_flowRate.dat');
a4 = load('bp1p7_flowRate.dat');
a5 = load('bp1p8_flowRate.dat');
a6 = load('bp1p9_flowRate.dat');

m = zeros(1,6);
m(1) = mean(a1(end-4:end,2));
m(2) = mean(a2(end-4:end,2));
m(3) = mean(a3(end-4:end,2));
m(4) = mean(a4(end-4:end,2));
m(5) = mean(a5(end-4:end,2));
m(6) = mean(a6(end-4:end,2));
p = linspace(1.4,1.9,6);
pp = p*1.4;
mm = m/(0.5*pi*0.1156);

%%
plot(pp,mm,'-or', 'Linewidth', 1.2);
grid on
%ylim([0.056, 0.065]);
ylim([0.0 0.35]);
xlim([2.0 3.4]);


ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
ax.TickLabelInterpreter = 'latex';
yticks([0.0 0.1 0.2 0.3]);
set(gcf, 'Position',  [100, 100, 800, 450]);
%title('Mass flow rate variation with time at three different locations for $p_b^* = 1.4$','Interpreter', 'latex');
xlabel('${p_b}/{p_i}$', 'Interpreter', 'latex');
ylabel('$\dot{m}_m$', 'Interpreter', 'latex');

%%
axes('Position',  [.4, .2, .3, .3]);
box on
plot(pp,mm,'-or', 'Linewidth', 1.2);
xlim([2.0 2.66]);
ylim([0.3 0.35]);
ax = gca;
ax.FontSize = 12;
ax.XAxis.LineWidth = 1.0;
ax.YAxis.LineWidth = 1.0;
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.TickLabelInterpreter = 'latex';