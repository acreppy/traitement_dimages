function [x, y, ux ,uy, num_app, tracking, pos_x, pos_y]=PTV_fluo_adama(root,i1,i2,inc)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dÃ©finition des paramÃ¨tres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_part_max=3000;
separation_max=10; %distance de sÃ©paration maximale autorisÃ©e

filename_result=[root '_results_adama.mat']; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gÃ©nÃ©ration des noms des fichiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp=0;
for i=i1:inc:i2
    cp=cp+1;
    listing_image(cp).name=[root num2str(i,'%04d') '.tif'];
end



x=zeros(cp-1,nb_part_max);
y=zeros(cp-1,nb_part_max);
ux=zeros(cp-1,nb_part_max);
uy=zeros(cp-1,nb_part_max);

num=zeros(cp,1);
tracking=zeros(cp-1,nb_part_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creation templace
% mask de taille (n,n) 
% gaussienne d'Ã©cart type sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=10;
sigma=3;
i=[1:n]-(n+1)/2;
[X Y]=meshgrid(i,i);
mask=exp(-(X.^2+Y.^2)/(2*sigma^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DÃ©tection des particules sur chacunes des images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:cp
%    name=[root num2str(i,'%04d') '.tif']
    name=listing_image(i).name;
    disp(name)
    a=double(imread(name));
    im=a(20:end,:);
    


    [x_temp y_temp num(i)]=detection(im,mask);
    pos_x(i,1:num(i))=x_temp;
    pos_y(i,1:num(i))=y_temp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appariement des particles entre deux images succesives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:cp-1
    num1=num(i);
    x1=pos_x(i,1:num1);
    y1=pos_y(i,1:num1);
    
    num2=num(i+1);
    x2=pos_x(i+1,1:num2);
    y2=pos_y(i+1,1:num2);
    
    
    % calcul de toutes les distances interparticulaires
    [X1 X2]=meshgrid(x1,x2);
    [Y1 Y2]=meshgrid(y1,y2);
    
    D=((X2-X1).^2+(Y2-Y1).^2);
    
    % recherche du plus proche voisin qui correspond au dÃ©placement minimum
    [dep1_2 part2]=min(D,[],1);
    dx=x2(part2)-x1;
    dy=y2(part2)-y1;
    d=(dx.^2+dy.^2).^0.5;
   
    
    % Selection des bons appariement 
    ind_vrai=find(d<separation_max); % critÃ¨re sur la sÃ©paration max 
    num_vrai=length(ind_vrai);
    x(i,1:num_vrai)=x1(ind_vrai);
    y(i,1:num_vrai)=y1(ind_vrai);
    ux(i,1:num_vrai)=dx(ind_vrai);
    uy(i,1:num_vrai)=dy(ind_vrai);
    tracking(i,ind_vrai)=part2(ind_vrai);
    num_app(i)=num_vrai;      
   
end

%x(i+1,1:num(i+1))=pos_x(i+1,1:num(i+1));
%y(i+1,1:num(i+1))=pos_y(i+1,1:num(i+1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sauvegarde des resulats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(filename_result,'x','y','ux','uy','num_app','tracking','pos_x', 'pos_y')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
