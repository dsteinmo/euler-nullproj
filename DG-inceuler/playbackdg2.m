%plotstart;

figure(1); clf;
%set(gcf,'Visible','Off');
set(gcf,'renderer','painters');

stuff = ls('out*.mat');
filelist = strread(stuff,'%s');
filelist = sort(filelist); %put into proper time ordering

numfiles = length(filelist);

for ii=1:numfiles
   filenamestr = filelist{ii};
   load(filenamestr);
   colormap(darkjet); pf2dquad(N,x,y,Qnp1(:,:,1));
   colorbar; 
   rho = Qnp1(:,:,1);
   axis equal; caxis([min(rho(:)) max(rho(:)) ]);
   title(['t=' num2str(time(end))]);
   drawnow; 
   ii
   %print('-dpng', ['frame' sprintf('%07d',ii) '.png']);
end
   
