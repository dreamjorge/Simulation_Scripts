%Script For Genereate Slices, we need run first MainHermite for obtain Wo
%Matrix
znorm = z;
xnorm = x;
ynorm = y;
xmin  = min(znorm(:));
xmax  = max(znorm(:)); 
ymin  = min(xnorm(:));  
ymax  = max(xnorm(:)); 
zmin  = min(xnorm(:));  
zmax  = max(xnorm(:));

% Introduce Nan Values on Marix
bb         = Wo;
iindex     = find(abs(bb)<=.0005);
bb(iindex) = NaN;


close(figure(101))
figure (101)
plotPropagatedRays(rayH11,rayH12,rayH21,rayH22,1,1)
hold on 
daspect([1.5,.1,.1]) 
axis tight 
view(44,16) 
camzoom(1.4) 
camproj perspective
axis off
%%Correction of axis on Slices
aa = permute(bb,[3 2 1]);
%distances of propagation
zprop  = [0, 6*Dz/3 , Dz/2, 0.99*Dz];
hzprop = slice(znorm,ynorm,xnorm,abs(aa ),zprop,[],[]); 
hzprop(1).FaceColor = 'interp'; 
hzprop(2).FaceColor = 'interp';
hzprop(3).FaceColor = 'interp';
hzprop(4).FaceColor = 'interp';
hzprop(1).EdgeColor = 'none'; 
hzprop(2).EdgeColor = 'none'; 
hzprop(3).EdgeColor = 'none'; 
hzprop(4).EdgeColor = 'none'; 
colormap(mapgreen)
hold off
export_fig('SelfHealingHermite3D','-png','-transparent')