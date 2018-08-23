 clear all
clc
format long

addpath C:\Users\will\Desktop\Misure\Misure-20180725

%%importare sweep per ottenere slope e f0

Sweep=importdata('20180725-1s6d5-15mW-150mVpk-sweepRBM3.txt',';',8);
figure(1)
title('s6d5 Optical modulation -RBM3-15mW-150mVpk','fontsize',30)
xAmpli=Sweep.data(1:end/2,1);
yAmpli=Sweep.data(1:end/2,2);
xPhase=Sweep.data(end/2+1:end,1);
yPhase=Sweep.data(end/2+1:end,2);
ax1 = subplot(2,1,1);
        hold on
        plot(xAmpli*1e-6,yAmpli,'linewidth',1.2 )
        ylabel(ax1,'Amplitude(dB)','fontsize',18)
        xlabel(ax1,'Frequency (MHz)','fontsize',18);
        
        box on
ax2 = subplot(2,1,2);
        hold on
        plot(xPhase*1e-6,yPhase,'linewidth',1.2 )  
        ylabel(ax2,'Phase(deg)','fontsize',18)
        xlabel(ax2,'Frequency (MHz)','fontsize',18);
        box on
[xP,yP]=ginput(3)
slope=(xP(3)-xP(2))/(yP(3)-yP(2));
f0=xP(1);

%%scrittura file con punti scelti , slope, f0
f = fopen('C:\Users\will\Desktop\Misure\Misure-20180725\s6d5 Optical modulation -RBM3-15mW-150mVpk.DAT','w')
fprintf(f,'%5s %f %f\n','P(1)',xP(1),yP(1));
fprintf(f,'%5s %f %f\n','P(2)',xP(2),yP(2));
fprintf(f,'%5s %f %f\n\n','P(3)',xP(3),yP(3));
fprintf(f,'%2s %f\n','f0',f0);
fprintf(f,'%5s %f\n','slope',slope);
fclose(f);

% importare file della traccia temporale della fase
PhaseTimeTrace=importdata('20180725-1s6d5-15mW-150mVpk-3Modes-openloop.txt',';',9);
Pos_fin=find(PhaseTimeTrace.data(:,1)==0); %%identifica la posizione nel vettore tempo di t=0 che corrisponde all'ultimo acquisito per ogni modo
n_modi=length(Pos_fin); %%identifica il numero di modi presenti nella traccia temporale in funzione del numero di zeri presenti
% % Values=importdata('s6d5 Optical modulation -RBM1-15mW-150mVpk.DAT');
% % f0=Values.data(4,1); %Hz
% % slope=Values.data(5,1); %Hz/deg
if slope<0
    slope=-slope;
end    

time=PhaseTimeTrace.data(1:Pos_fin(1),1);
N=numel(time);
for m=1:n_modi 
    clear phase
    phase(:)=PhaseTimeTrace.data(1+(m-1)*N:Pos_fin(m),2);
    tau=zeros(floor(N/2),1);
    p=0;

    for n=1:floor(N/2)
        tau(n)=n*(time(2)-time(1));
        s=1;
        diff=0;
        for k=1:floor(N/n)
            pmean=mean(phase((s):(n+s-1)));
            diff=diff+(pmean-p)^2;
            p=pmean;
            s=s+n;
        end
        Y(n)=sqrt(1/(2*(floor(N/n)-1))*diff);
    end
    Allan(:,m)=Y;
    Allan(:,m)=Y*slope/f0;
end
%scrittura file con x=tempo di integrazione e y=Allan deviation
dlmwrite('C:\Users\will\Desktop\Misure\Misure-20180725\Allan-20180725-15mW-150mVpk-OpenLoop.DAT',[tau Allan])

figure(2)
loglog(tau,Allan,'o');
xlabel('time (sec)');
ylabel('\sigma_{Allan}(\tau)')
set(gca,'fontsize',16)
title('20180725-15mW-150mVpk-OpenLoop','fontsize',12)
% legend({'slope=' num2str(slope)},'FontSize',16);
%legend( sprintf('slope %f Hz/deg', slope) )

box on
%print('C:\Users\will\Desktop\Misure\Misure-20180718\RBM1\20180718-15mW-100mVpk-sweep-1526~5nm-RBM1SecondPeak-1~717KHz-1min-Figure','-dpng')
savefig('C:\Users\will\Desktop\Misure\Misure-20180718\RBM1\20180718-15mW-150mVpk-sweep-1526~5nm-RBM1-3~43KHz-1min-Figure.fig')