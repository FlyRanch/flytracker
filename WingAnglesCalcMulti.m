close all; clear all; clc

expname = ['WingAngles_exp098c';
           'WingAngles_exp101c';
           'WingAngles_exp099c';
           'WingAngles_exp100c';
           'WingAngles_exp102c';
           'WingAngles_exp104c';
           'WingAngles_exp083c';
           'WingAngles_exp035c';
           ];
expnum  = [36;
           39;
           37;
           38;
           40;
           42;
           24;
            6;
          ];
wbeats  = [ 5 10  14 17 22 26 31  0  0  0;
            5 10  14 17 21 25 28 32 35 40;
            3  7  10 14 18 21 25 29 33 36;
            6 10  14 18  0  0  0  0  0  0;
           20 24  27 31 35 39 43 46 49 53;
           22 25  30 34  0  0  0  0  0  0; 
           1   4   6  8 10 14  0  0  0  0;
           ]; 
       
stp = 1/6;
wbeats = wbeats./stp;
wbeatsnum = [ 7;
             10;
             10;
              4;
             10;
              4;
              6;
              ];
color = ['r:';'b:';'g:';'m:';'c:';'r:';'b:';'g:'];
          
          
load gmc/Combo2Results_n.mat
      
scrsz = get(0,'ScreenSize');

SAmp = figure('Position',[1 scrsz(4)*0.6 scrsz(3)/6 scrsz(4)*0.4]); 
SPD  = figure('Position',[1 1 scrsz(3)/6 scrsz(4)*0.4]);
BO   = figure('Position',[scrsz(3)/6 1  scrsz(3)/6 scrsz(4)*0.4]);     
AoA  = figure('Position',[scrsz(3)/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.4]); 
WBF  = figure('Position',[scrsz(3)*2/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.4]); 
Amp  = figure('Position',[scrsz(3)*2/6 1  scrsz(3)/6 scrsz(4)*0.4]);
%GMC  = figure('Position',[scrsz(3)*3/6 scrsz(4)*0.6  scrsz(3)/6 scrsz(4)*0.4]); 

for ij = 1

    load(expname(ij,:))

%%  %===================
    %Plot wing angles
    %t_end = movidx(end);

    %===================
    %Plot Stroke Amplitude   

    %tt = movidx;
    t_end = round(tt(end));
    figure(SAmp); 
    plot(tt,phi_R*(180/pi),'r-','linewidth',2);
    hold on
    plot(tt,phi_L*(180/pi),'b-','linewidth',2);
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-110 110],'ytick',-110:55:110,'yticklabel',{'-110','','0','','110'});
    hold off

    %===================
    %Plot Stroke Plane Deviation (w.r.t. 62deg)    
    figure(SPD); 
    plot(tt,theta_R*(180/pi),'r-','linewidth',2);
    hold on
    plot(tt,theta_L*(180/pi),'b-','linewidth',2);
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-40 40],'ytick',-40:20:40,'yticklabel',{'-40','','0','','40'});
    hold off
    
    %===================
    %Plot Angle of Attack (geometric)    
	figure(AoA)
    plot(tt(1:end-1),alpha_R(1:end-1)*(180/pi),'r-','linewidth',2);
    hold on
    plot(tt(1:end-1),alpha_L(1:end-1)*(180/pi),'b-','linewidth',2);
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    set(gca,'ylim',[-180 180],'ytick',-180:90:180,'yticklabel',{'-180','','0','','180'});
    hold off
    
    %===================
    %Plot Body orientation
    figure(BO);
    plot(tt,BodyAng_auto(1,:)*(180/pi),'r-','linewidth',2);
    hold on;
    plot(tt,BodyAng_auto(2,:)*(180/pi),'g-','linewidth',2);
    plot(tt,(BodyAng_auto(3,:) - BodyAng_auto(3,1))*(180/pi),'b-','linewidth',2);
    set(gca,'plotboxaspectratio',[3 1 1],'xlim',[tt(1) t_end]);
    box off
    hold off

    %===================
    %Wingbeat Frequency
    [maxAmpR, minAmpR] = peakdet(phi_R, 0.1, tt);
    [maxAmpL, minAmpL] = peakdet(phi_L, 0.1, tt);
   
    figure(SAmp); hold on
    plot(minAmpR(:,1), minAmpR(:,2).*180/pi, 'g*');
    plot(minAmpL(:,1), minAmpL(:,2).*180/pi, 'y*');  
    plot(maxAmpR(:,1), maxAmpR(:,2).*180/pi, 'r*');
    plot(maxAmpL(:,1), maxAmpL(:,2).*180/pi, 'm*');
    hold off
    
    for ik=1:length(minAmpR(:,2))-1
        wbfR(1,ik) = 1/(minAmpR(ik+1,1)-minAmpR(ik,1))*1000;
    end

     for ik=1:length(minAmpL(:,2))-1
        wbfL(1,ik) = 1/(minAmpL(ik+1,1)-minAmpL(ik,1))*1000;
     end
    
    wbfR = [wbfR(1,1) wbfR];
    wbfL = [wbfL(1,1) wbfL];
     
    figure(WBF); 
    plot(minAmpR(:,1)-minAmpR(1,1),wbfR,color(ij,:),'LineWidth',2); hold on   
    plot(minAmpR(:,1)-minAmpR(1,1),wbfR,'kx','LineWidth',2)
 
    plot(minAmpL(:,1)-minAmpR(1,1),wbfL,color(ij,:),'LineWidth',2)
    plot(minAmpL(:,1)-minAmpR(1,1),wbfL,'ko','LineWidth',1.2)    
    set(gca,'plotboxaspectratio',[3 1 1]);
    set(gca,'xlim',[0 40])
    %set(gca,'ylim',[0 300],'ytick',0:100:300,'yticklabel',{'0','','','300'});
    
    %if they occurr subsequently in time
    for iii = 1:length(minAmpR(:,1))-1
        for jjj = 1:length(maxAmpR(:,1))
            if (maxAmpR(jjj,1)>minAmpR(iii,1))&&(maxAmpR(jjj,1)<minAmpR(iii+1,1))
                A(iii) = (maxAmpR(jjj,2) - minAmpR(iii,2)).*180/pi;
                break;
            end
        end
    end 

    for iii = 1:length(minAmpL(:,1))-1
        for jjj = 1:length(maxAmpL(:,1))
            if (maxAmpL(jjj,1)>minAmpL(iii,1))&&(maxAmpL(jjj,1)<minAmpL(iii+1,1))
                B(iii) = (maxAmpL(jjj,2) - minAmpL(iii,2)).*180/pi;
                break;
            end
        end
    end 
    
    figure(Amp)
    plot(1:length(A),A,color(ij,:),'LineWidth',2);hold on
    plot(1:length(A),A,'kx','LineWidth',2)
    plot(1:length(B),B,color(ij,:),'LineWidth',2)
    plot(1:length(B),B,'ko','LineWidth',1.2)    
    set(gca,'plotboxaspectratio',[3 1 1]);
    set(gca,'xtick',1:1:length(A),'xticklabel',{'1','2','3','4','5','6','7','8','9'});
    box off
    
    clear wbfR wbfL
   
%    idbeg = find(allData(expnum(ij)).f_align==0);
%    idend = length(allData(expnum(ij)).f_align);
    

    %figure(GMC);
%    figure(BO); hold on
    %gtime = allData(expnum(ij)).f_align(idbeg:idend)./6;
    %plot(gtime,allData(expnum(ij)).exp.rpy(1,idbeg:idend).*-180/pi,'r:','linewidth',2);
    %hold on
    %plot(gtime,allData(expnum(ij)).exp.rpy(2,idbeg:idend).*180/pi,'g:','linewidth',2);
    %plot(gtime,allData(expnum(ij)).exp.rpy(3,idbeg:idend).*180/pi,'b:','linewidth',2);  
    %hold off
    %set(gca,'plotboxaspectratio',[3 1 1]);

end

figure(SAmp); 
title('S.Amp');box off  

figure(SPD);
title('S.P.D');box off

figure(BO);     
title('B.O:E - roll(r), pitch(g), yaw(b)');box off

figure(AoA); 
title('A-o-A');box off

figure(WBF);
title('WBF');box off

%figure(GMC);
%title('B.O:G - idem');box off