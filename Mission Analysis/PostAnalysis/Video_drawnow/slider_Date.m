% x = [-1 , 1/3, 1/3, 1, 1/3, 1/3,-1 ];
% y = [-1/3,-1/3,-1/2, 0, 1/2, 1/3, 1/3];
x = [-1/3, -1/3, -2/3, 0, 2/3, 1/3, 1/3 ];
y = [1.5, 0.5, 0.5, 0, 0.5, 0.5, 1.5];
g = hgtransform;
patch('XData',x,'YData',y,'FaceColor','yellow','Parent',g)

cs1 = '\begin{tabular}{lllllll} 2027 & 2028 & 2029 & 2030 & 2031 & 2032 & 2033 \end{tabular}';

axis equal
xlim([-10 10])
ylim([-10 10])

% By default, the units are normalized to the figure. The lower left corner of the figure maps to (0,0) and the upper right corner maps to (1,1)
annotation('rectangle',[.3 .8 .5 .001])
    
h_textbox1 = annotation('textbox', [0.25, 0.7, 0.5, 0.1],'Fontsize',11,'Interpreter', 'latex','EdgeColor','none');
pt1 = [-5 8 0];
pt2 = [5 8 0];
for t=linspace(0,1,100)
  g.Matrix = makehgtform('translate',pt1 + t*(pt2-pt1));
  hold on 
  set(h_textbox1, 'String', cs1);
  drawnow
end

%%
s1 = 1/2;
s2 = 2;
r1 = 0;
r2 = 2*pi/3;
for t=linspace(0,1,100)
  g.Matrix = makehgtform('translate',pt1 + t*(pt2-pt1), ...
                         'scale',s1 + t*(s2-s1), ...
                         'zrotate',r1 + t*(r2-r1));
  drawnow
end