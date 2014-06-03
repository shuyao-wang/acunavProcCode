%%

one = load('/emfd/vrk4/AcuNav/P1280104/hARFI_SWEI/res_20140530135359');
two = load('/emfd/vrk4/AcuNav/P1280104/hARFI_SWEI/res_20140530135416');

figure
set(gcf,'Position',[1630 1028 224 582])
for i=10:40
    for j=1:41
        subplot(121)
        imagesc([],one.radial,one.arfidata(:,:,i,j),[0 4])
        title(['Set 1: ',num2str(i),' ',num2str(j)])
        subplot(122)
        imagesc([],two.radial,two.arfidata(:,:,i,j),[0 4])
        title(['Set 2: ',num2str(i),' ',num2str(j)])
        pause
    end
end
