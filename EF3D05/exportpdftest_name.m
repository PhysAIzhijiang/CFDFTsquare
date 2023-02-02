function []=exportpdftest_name(input_text)
h = gcf;
%plot(1:10);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,input_text,'-dpdf','-r0')