close all; clear all; clc

expname = ['WingAngles_exp098c';
           'WingAngles_exp101c';
           'WingAngles_exp099c';
           'WingAngles_exp102c';
           'WingAngles_exp104c';           
           'WingAngles_exp095c';           
           %'WingAngles_exp083c';
           'WingAngles_exp035c';
           'WingAngles_exp100c';    %this is the shortest one, so exclude plotting tail           
           ];
expnum  = [36;
           39;
           37;
           40;
           42;           
           33;           
           %24;
            6;
           38;  %this is the shortest one, so exclude plotting tail            
          ];
nofexp = 8;
color = ['r--';'b--';'g--';'m--';'c--';'r--';'b--';'g--'];
if 0
    expname = ['WingAngles_exp083c';
               'WingAngles_exp035c';
               ];
    expnum  = [ 24;
                6;
              ];
    nofexp = 2;
end      
load gmc/Combo2Results_n.mat
      
scrsz = get(0,'ScreenSize');

test = figure;

%MHD  = figure('Position',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.3 scrsz(4)*0.6]);

weightvector = zeros(1,1000);
sumvector = zeros(1,1000);


%%

for ij = 1:nofexp
    load(expname(ij,:))

    for icr = 1:length(BodyAng_auto(1,:))-1
        tmpBA(1,icr) = (BodyAng_auto(1,icr+1) - BodyAng_auto(1,icr))/(1/6);
        tmpBA(2,icr) = (BodyAng_auto(2,icr+1) - BodyAng_auto(2,icr))/(1/6);
        tmpBA(3,icr) = (BodyAng_auto(3,icr+1) - BodyAng_auto(3,icr))/(1/6);

    end
    rBodyAng_auto(1,:) = [ 0 tmpBA(1,:) ] ; 
    rBodyAng_auto(2,:) = [ 0 tmpBA(2,:) ] ;       
    rBodyAng_auto(3,:) = [ 0 tmpBA(3,:) ] ;                
        
    flystorer(ij) = {[tt;                    %01
                      phi_R;                 %02
                      phi_L;                 %03
                      theta_R;               %04
                      theta_L;               %05
                      alpha_R;               %06
                      alpha_L;               %07
                      rBodyAng_auto(1,:);    %08
                      rBodyAng_auto(2,:);    %09
                      rBodyAng_auto(3,:);    %10
                      BodyAng_auto(1,:);     %11
                      BodyAng_auto(2,:);     %12
                      BodyAng_auto(3,:)]};   %13                  
                  
    clear rBodyAng_auto tmpBA
    %===================
    %Get the peaks
    [maxAmpR, minAmpR] = peakdet(phi_R, 0.5, tt);
    [maxAmpL, minAmpL] = peakdet(phi_L, 0.5, tt);

    %Get indices
    indicesR{ij,:} = round(minAmpR(:,1).*6);
    indicesL{ij,:} = round(minAmpL(:,1).*6);
    
%    splineMe{ij,:} = alpha_L(indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe{ij,:} = flystorer{ij}(2,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe2{ij,:} = flystorer{ij}(3,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe3{ij,:} = flystorer{ij}(4,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe4{ij,:} = flystorer{ij}(5,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe5{ij,:} = flystorer{ij}(6,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe6{ij,:} = flystorer{ij}(7,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);    
    splineMe7{ij,:} = flystorer{ij}(8,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    
    
    
    splineMe8{ij,:} = flystorer{ij}(9,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);
    splineMe9{ij,:} = flystorer{ij}(10,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);              
    splineMe10{ij,:} = flystorer{ij}(11,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);   
    splineMe11{ij,:} = flystorer{ij}(12,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);   
    splineMe12{ij,:} = flystorer{ij}(13,indicesR{ij,1}(1,1)+1:indicesR{ij,1}(4,1)+1);       
    tmp(ij,:) = (indicesR{ij,1}(4,1) - indicesR{ij,1}(1,1))+1;
    
    figure(test); 
    plot(BodyAng_auto(1,:).*(180/pi),'r');hold on;
    plot(BodyAng_auto(2,:).*(180/pi),'b');
    plot(BodyAng_auto(3,:).*(180/pi),'g');
end

spln2 = max(tmp);

for ik = 1:nofexp
    spln1(ik,:) = length(splineMe{ik,1});
    splineD{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD2{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe2{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD3{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe3{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD4{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe4{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD5{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe5{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD6{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe6{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD7{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe7{ik,1}(1,:),1:1/spln2:2-1/spln2);    
    splineD8{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe8{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD9{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe9{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD10{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe10{ik,1}(1,:),1:1/spln2:2-1/spln2);    
    splineD11{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe11{ik,1}(1,:),1:1/spln2:2-1/spln2);
    splineD12{ik,:} = spline(1:1/spln1(ik,:):2-1/spln1(ik,:),splineMe12{ik,1}(1,:),1:1/spln2:2-1/spln2);    
end

for ik = 1:nofexp
    splineM(ik,:)  = splineD{ik,:}(1,:);
    splineM2(ik,:) = splineD2{ik,:}(1,:);
    splineM3(ik,:) = splineD3{ik,:}(1,:);
    splineM4(ik,:) = splineD4{ik,:}(1,:);
    splineM5(ik,:) = splineD5{ik,:}(1,:);
    splineM6(ik,:) = splineD6{ik,:}(1,:);
    splineM7(ik,:) = splineD7{ik,:}(1,:);
    splineM8(ik,:) = splineD8{ik,:}(1,:);
    splineM9(ik,:) = splineD9{ik,:}(1,:);    
    splineM10(ik,:) = splineD10{ik,:}(1,:);
    splineM11(ik,:) = splineD11{ik,:}(1,:);
    splineM12(ik,:) = splineD12{ik,:}(1,:);     
end


figure('Position',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.3 scrsz(4)*0.8]);

subplot(6,1,1);plot(1:3/(spln2-1):4,splineM.*180/pi,'b'); hold on; title('Right & Rates')
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM).*180/pi,'r','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM).*180/pi + std(splineM).*180/pi,'r--','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM).*180/pi - std(splineM).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,1);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(2,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,1);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(2,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,2);plot(1:3/(spln2-1):4,splineM3.*180/pi,'b'); hold on
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM3).*180/pi,'r','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM3).*180/pi + std(splineM3).*180/pi,'r--','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM3).*180/pi - std(splineM3).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,2);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(4,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,2);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(4,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,3);plot(1:3/(spln2-1):4,splineM5.*180/pi,'b'); hold on
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM5).*180/pi,'r','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM5).*180/pi + std(splineM5).*180/pi,'r--','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM5).*180/pi - std(splineM5).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,3);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(6,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,3);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(6,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,4);plot(1:3/(spln2-1):4,splineM7.*180/pi,'b'); hold on
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM7).*180/pi,'r','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM7).*180/pi + std(splineM7).*180/pi,'r--','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM7).*180/pi - std(splineM7).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,4);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(8,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,4);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(8,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,5);plot(1:3/(spln2-1):4,splineM8.*180/pi,'b'); hold on
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM8).*180/pi,'r','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM8).*180/pi + std(splineM8).*180/pi,'r--','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM8).*180/pi - std(splineM8).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,5);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(9,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,5);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(9,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,6);plot(1:3/(spln2-1):4,splineM9.*180/pi,'b'); hold on
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM9).*180/pi,'r','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM9).*180/pi + std(splineM9).*180/pi,'r--','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM9).*180/pi - std(splineM9).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,6);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(10,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,6);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(10,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;


%Left Wing

figure('Position',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.3 scrsz(4)*0.8]);

subplot(6,1,1);plot(1:3/(spln2-1):4,splineM2.*180/pi,'b'); hold on; title('Left & Rates')
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM2).*180/pi,'r','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM2).*180/pi + std(splineM2).*180/pi,'r--','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM2).*180/pi - std(splineM2).*180/pi,'r--','linewidth',2)
set(gca,'xtick',[],'xticklabel','');box off;
%
for plotik = 1:nofexp
subplot(6,1,1);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(3,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,1);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(3,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,2);plot(1:3/(spln2-1):4,splineM4.*180/pi,'b'); hold on
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM4).*180/pi,'r','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM4).*180/pi + std(splineM4).*180/pi,'r--','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM4).*180/pi - std(splineM4).*180/pi,'r--','linewidth',2)
set(gca,'xtick',[],'xticklabel','');box off;
%
for plotik = 1:nofexp
subplot(6,1,2);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(5,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,2);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(5,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,3);plot(1:3/(spln2-1):4,splineM6.*180/pi,'b'); hold on
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM6).*180/pi,'r','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM6).*180/pi + std(splineM2).*180/pi,'r--','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM6).*180/pi - std(splineM2).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,3);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(7,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,3);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(7,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,4);plot(1:3/(spln2-1):4,splineM7.*180/pi,'b'); hold on
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM7).*180/pi,'r','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM7).*180/pi + std(splineM7).*180/pi,'r--','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM7).*180/pi - std(splineM7).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,4);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(8,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,4);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(8,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,5);plot(1:3/(spln2-1):4,splineM8.*180/pi,'b'); hold on
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM8).*180/pi,'r','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM8).*180/pi + std(splineM8).*180/pi,'r--','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM8).*180/pi - std(splineM8).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,5);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(9,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,5);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(9,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,6);plot(1:3/(spln2-1):4,splineM9.*180/pi,'b'); hold on
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM9).*180/pi,'r','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM9).*180/pi + std(splineM9).*180/pi,'r--','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM9).*180/pi - std(splineM9).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,6);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(10,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,6);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(10,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

%Rite wing mate

figure('Position',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.3 scrsz(4)*0.8]);

subplot(6,1,1);plot(1:3/(spln2-1):4,splineM.*180/pi,'b'); hold on; title('Right & Angles')
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM).*180/pi,'r','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM).*180/pi + std(splineM).*180/pi,'r--','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM).*180/pi - std(splineM).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,1);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(2,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp-1
subplot(6,1,1);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(2,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,2);plot(1:3/(spln2-1):4,splineM3.*180/pi,'b'); hold on
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM3).*180/pi,'r','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM3).*180/pi + std(splineM3).*180/pi,'r--','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM3).*180/pi - std(splineM3).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,2);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(4,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,2);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(4,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,3);plot(1:3/(spln2-1):4,splineM5.*180/pi,'b'); hold on
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM5).*180/pi,'r','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM5).*180/pi + std(splineM5).*180/pi,'r--','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM5).*180/pi - std(splineM5).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,3);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(6,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1 
subplot(6,1,3);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(6,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,4);plot(1:3/(spln2-1):4,splineM10.*180/pi,'b'); hold on
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM10).*180/pi,'r','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM10).*180/pi + std(splineM10).*180/pi,'r--','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM10).*180/pi - std(splineM10).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,4);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(11,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,4);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(11,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,5);plot(1:3/(spln2-1):4,splineM11.*180/pi,'b'); hold on
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM11).*180/pi,'r','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM11).*180/pi + std(splineM11).*180/pi,'r--','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM11).*180/pi - std(splineM11).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,5);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(12,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,5);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(12,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,6);plot(1:3/(spln2-1):4,splineM12.*180/pi,'b'); hold on
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM12).*180/pi,'r','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM12).*180/pi + std(splineM12).*180/pi,'r--','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM12).*180/pi - std(splineM12).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,6);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(13,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,6);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(13,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;


%Left Wing

figure('Position',[scrsz(3)*0.1 scrsz(4)*0.1 scrsz(3)*0.3 scrsz(4)*0.8]);

subplot(6,1,1);plot(1:3/(spln2-1):4,splineM2.*180/pi,'b'); hold on; title('Left & Angles')
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM2).*180/pi,'r','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM2).*180/pi + std(splineM2).*180/pi,'r--','linewidth',2)
subplot(6,1,1);plot(1:3/(spln2-1):4,mean(splineM2).*180/pi - std(splineM2).*180/pi,'r--','linewidth',2)
set(gca,'xtick',[],'xticklabel','');box off;
%
for plotik = 1:nofexp
subplot(6,1,1);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(3,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,1);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(3,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,2);plot(1:3/(spln2-1):4,splineM4.*180/pi,'b'); hold on
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM4).*180/pi,'r','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM4).*180/pi + std(splineM4).*180/pi,'r--','linewidth',2)
subplot(6,1,2);plot(1:3/(spln2-1):4,mean(splineM4).*180/pi - std(splineM4).*180/pi,'r--','linewidth',2)
set(gca,'xtick',[],'xticklabel','');box off;
%
for plotik = 1:nofexp
subplot(6,1,2);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(5,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,2);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(5,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,3);plot(1:3/(spln2-1):4,splineM6.*180/pi,'b'); hold on
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM6).*180/pi,'r','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM6).*180/pi + std(splineM2).*180/pi,'r--','linewidth',2)
subplot(6,1,3);plot(1:3/(spln2-1):4,mean(splineM6).*180/pi - std(splineM2).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,3);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(7,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,3);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(7,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,4);plot(1:3/(spln2-1):4,splineM10.*180/pi,'b'); hold on
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM10).*180/pi,'r','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM10).*180/pi + std(splineM10).*180/pi,'r--','linewidth',2)
subplot(6,1,4);plot(1:3/(spln2-1):4,mean(splineM10).*180/pi - std(splineM10).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,4);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(11,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,4);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(11,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,5);plot(1:3/(spln2-1):4,splineM11.*180/pi,'b'); hold on
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM11).*180/pi,'r','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM11).*180/pi + std(splineM11).*180/pi,'r--','linewidth',2)
subplot(6,1,5);plot(1:3/(spln2-1):4,mean(splineM11).*180/pi - std(splineM11).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,5);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(12,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,5);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(12,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;

subplot(6,1,6);plot(1:3/(spln2-1):4,splineM12.*180/pi,'b'); hold on
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM12).*180/pi,'r','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM12).*180/pi + std(splineM12).*180/pi,'r--','linewidth',2)
subplot(6,1,6);plot(1:3/(spln2-1):4,mean(splineM12).*180/pi - std(splineM12).*180/pi,'r--','linewidth',2)
%
for plotik = 1:nofexp
subplot(6,1,6);plot( 0:1/indicesR{plotik,1}(1,1):1,...
    flystorer{1,plotik}(13,( 1:indicesR{plotik,1}(1,1)+1 )).*180/pi,'k-');
end

for plotik = 1:nofexp - 1
subplot(6,1,6);plot( 4:1/( indicesR{plotik,1}(5,1)-indicesR{plotik,1}(4,1)-1 ):5,...
    flystorer{1,plotik}(13,( indicesR{plotik,1}(4,1)+1:indicesR{plotik,1}(5,1) )).*180/pi,'k-');
end
%
set(gca,'xtick',[0 1 2 3 4],'xticklabel',{'0' '1' '2' '3' '4'});box off;