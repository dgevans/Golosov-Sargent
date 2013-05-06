x = linspace(1,15,128);
y = linspace(1,15,128);
figure;
contour(x,y,flipud(P));
axis square;
colormap('Copper');
xlabel('X Distance(in mm)');
ylabel('Y Distance(in mm)'