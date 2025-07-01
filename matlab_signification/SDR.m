%% SDR simplex 

%% SDR simplex calculation
% 
% 
% <<SDR_kova.jpg>>
% 
% 
% In general, any data matrix containing presence-absence or abundance data can
% be analyzed.
% The function accepts species as row and sites as
% columns respectively the matrix shoudl be named as :data.
% 
%% Backgound Notes
% A conceptual framework is proposed to evaluate the relative importance of beta diversity, nestedness and agreement in species
% richness in presence?absence data matrices via partitioning pairwise gamma diversity into additive components. This
% is achieved by calculating three complementary indices that measure similarity, relative species replacement, and relative
% richness difference for all pairs of sites, and by displaying the results in a two-dimensional simplex diagram, or ternary plot.
% By summing two terms at a time, three one-dimensional simplices are derived correspondig to different contrasts: beta
% diversity versus similarity, species replacement versus nestedness and, finally, richness difference versus richness agreement.
% The simplex diagrams can be used to interpret underlying data structures by showing departure from randomness towards
% % well-interpretable directions, as demonstrated by artificial and actual examples. In particular, one may appreciate how
% far data structure deviates from three extreme model situations: perfect nestedness, anti-nestedness and perfect gradient.
% Throughout the paper, we pay special attention to the measurement and interpetation of beta diversity and nestedness for
% pairs of sites, because these concepts have been in focus of ecological reseach for decades. The novel method can be used in
% community ecology, conservation biology, and biogeography, whenever the objective is to recover explanatory ecological
% processes behind patterns conveyed by presence–absence data.
% 
%
%% Transform the matrix into binary data

[sor,oszlop]=size(data);
%%
bin=zeros(sor,oszlop);
for j=1:oszlop
    for i=1:sor
        if data(i,j)>0 
            bin(i,j)=1;
        else bin(i,j)=0;
        end
    end
end

%% Calculation of the SDR indexes

a=zeros((oszlop*(oszlop-1)/2),9);
k=1;

for i=1:oszlop
    
    for j=i+1:oszlop
        
        a(k,1)=sum((bin(:,i)+bin(:,j))==2);              %a
        a(k,2)=sum(bin(:,i)-bin(:,j)==1);               %b
        a(k,3)=sum(bin(:,i)-bin(:,j)==-1);                %c
        a(k,4)=sum((bin(:,i)+bin(:,j))~=0);               %n
        a(k,5)=a(k,1)/a(k,4);                    %S=a/n
        a(k,6)=abs(a(k,2)-a(k,3))/a(k,4);  % D=abs(b-c)/n
        a(k,7)=2*min(a(k,2),a(k,3))/a(k,4); %R=2min(b,c)/n
        a(k,8)=1/2*(2*a(k,5)+a(k,7))/(a(k,5)+a(k,6)+a(k,7)) ; % x 
        a(k,9)=sqrt(3)/2*a(k,7)/(a(k,5)+a(k,6)+a(k,7)) ;   %y
        k=k+1;
    end
end
%% Triangular figure

harom=[0,0;1,0;0.500000000000000,0.866025403784439];

 fill(harom(:,1),harom(:,2),'w')
 hold
  line([0.25 0.5],[0.43 0.33])
  line([0.75 0.5],[0.43 0.33])
  line([0.5 0.5],[00 0.33])
 scatter(a(:,8),a(:,9),'filled')
 hold off
 %% SDR vales
 % The output sequence by columns is R, D and S

 sdr(:,1)=mean(a(:,5));
 sdr(:,2)=mean(a(:,6));
 sdr(:,3)=mean(a(:,7));
 sdr
 %% References
 % Podani J. and D. Schmera. 2011. A new conceptual and methodological framework for exploring and explaining pattern in presence-absence data. – Oikos, 120, 1625–1638
 % Oikos 120: 1625–1638, 2011
%  doi: 10.1111/j.1600-0706.2011.19451.x.
