% Detection des  des particules :  correlation avec mask 
% Localisation des particules : regresssion gaussienne
%

function [x y,num]=detection(im,mask)
ni_mask=size(mask,1);
nj_mask=size(mask,2);

% evaluation de la matrice  de correlation et seuillage
c2=normxcorr2(mask,im);
nic2=size(c2,1);
njc2=size(c2,2);
b2=im2bw(c2,0.7);

% [i2 j2]=find(b2);
% i2=i2-(n-1)/2;
% j2=j2-(n-1)/2;
%hold on; plot(j2,i2,'og')




% labelisation des particules
CC = bwconncomp(b2);
num=CC.NumObjects;
L = labelmatrix(CC);




% extraction centre de masse pour regression gaussienne

STATS = regionprops(b2,c2,'WeightedCentroid');  
centroids = cat(1, STATS.WeightedCentroid);



% DÃ©termination de la position des particules par 
% Regression gaussienne effectuÃ©ee sur la matrice de correlation

for part=1:num
    
% recherche de la position en  pixel du maximum de chacune des particules dÃ©tectÃ©e 
%     ind=find(L==part) ;
%     i0=find(c2(ind)==max(c2(ind)));
%     ind_max=ind(i0);
%     j_max= ceil(ind_max/nic2);
%     i_max= ind_max-nic2*(j_max-1);

    i_max=round(centroids(part,2));
    j_max=round(centroids(part,1));

    % extraction de la position par rÃ©gression gaussienne

    % Nobach, H., & Honkanen, M. (2005).
    % Two-dimensional Gaussian regression for sub-pixel displacement estimation
    % in particle image velocimetry or particle position estimation in particle tracking velocimetry.
    % Experiments in Fluids, 38(4), 511â€“515.
    c10=0;
    c01=0;
    c11=0;
    c20=0;
    c02=0;
    c00=0;
    for i=-1:1
        for j=-1:1
            temp=log(c2(i_max+i,j_max+j));
            c10=c10+i*temp;
            c01=c01+j*temp;
            c11=c11+i*j*temp;
            c20=c20+(3*i^2-2)*temp;
            c02=c02+(3*j^2-2)*temp;
            c00=c00+(5-3*i^2-3*j^2)*temp;
        end
    end
    c10=c10/6;
    c01=c01/6;
    c11=c11/4;
    c20=c20/6;
    c02=c02/6;
    c00=c00/9;

    dx=(c11*c01-2*c10*c02)/(4*c20*c02-c11^2);
    dy=(c11*c10-2*c01*c20)/(4*c20*c02-c11^2);

    center_i(part)=i_max+dx;
    center_j(part)=j_max+dy;
    
end
x=center_j(:)-(nj_mask-1)/2; 
y=center_i(:)-(ni_mask-1)/2;