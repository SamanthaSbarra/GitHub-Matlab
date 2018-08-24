clear all
clc
format long

%%introdurre la cartella di lavoro
percorso='C:\Users\will\Dropbox\Misure\Misure-20180725';
cd 'C:\Users\will\Dropbox\Misure\Misure-20180725'
fileTxt = dir('*.txt');
n_modi=0;
for i=1:length(fileTxt)
    k=strfind(fileTxt(i).name,'sweep');
    if ~isempty(k)
        n_modi=n_modi+1;
        k=strfind(fileTxt(i).name,'RBM');
        j=str2num(fileTxt(i).name(k+1));
        Sweep(n_modi)=importdata(fileTxt(i).name,';');
    end
end    

for i=1:length(fileTxt)
    k=strfind(fileTxt(i).name,'PLL');
    z=strfind(fileTxt(i).name,'OpenLoop');
    if ~isempty(k)||~isempty(z)
        TracciaTemporale=importdata(fileTxt(i).name,';');
        s=strfind(fileTxt(i).name,'s');
        for w=1:length(s)
            if fileTxt(i).name(s(w)+2)=='d'
                serie=fileTxt(i).name(s(w)+1);
                disco=fileTxt(i).name(s(w)+3);
            end
        end    
        power=strfind(fileTxt(i).name,'mW');
        power=fileTxt(i).name(power-2:power+1);
        amplit=strfind(fileTxt(i).name,'mVpk');
        amplit=fileTxt(i).name(amplit-3:amplit+3);
        date=strfind(fileTxt(i).name,'2018');
        date=fileTxt(i).name(date:date+7);
        if ~isempty(z)>0
            Mode='OpenLoop';
        else    
            Mode='PLL';
        end
    end   
end    

for i=1:n_modi
    figure(i)
    xAmpli=Sweep(i).data(1:end/2,1);
    yAmpli=Sweep(i).data(1:end/2,2);
    xPhase=Sweep(i).data(end/2+1:end,1);
    yPhase=Sweep(i).data(end/2+1:end,2);
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
    slope(i)=(xP(3)-xP(2))/(yP(3)-yP(2));
    f0(i)=xP(1);
    close
    %%scrittura file con punti scelti , slope, f0. Inserire Percorso di
    %%scrittura
    fileName=sprintf('s%sd%s-RBM%d-%s-%s.DAT',serie,disco,i,power,amplit);
    f = fopen(percorso+"\"+fileName,'w');
    fprintf(f,'%5s %f %f\n','P(1)',xP(1),yP(1));
    fprintf(f,'%5s %f %f\n','P(2)',xP(2),yP(2));
    fprintf(f,'%5s %f %f\n\n','P(3)',xP(3),yP(3));
    fprintf(f,'%2s %f\n','f0',f0(i));
    fprintf(f,'%5s %f\n','slope',slope(i));
    fclose(f);
end
% importare file della traccia temporale della fase
Pos_fin=find(TracciaTemporale.data(:,1)==0); %%identifica la posizione nel vettore tempo di t=0 che corrisponde all'ultimo acquisito per ogni modo
n_modi=length(Pos_fin); %%identifica il numero di modi presenti nella traccia temporale in funzione del numero di zeri presenti   

time=TracciaTemporale.data(1:Pos_fin(1),1);
N=numel(time);
for m=1:n_modi 
    clear phase
    phase(:)=TracciaTemporale.data(1+(m-1)*N:Pos_fin(m),2);
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
    if slope(m)<0
        slope(m)=-slope(m);
    end 
    Allan(:,m)=Y*slope(m)/f0(m);
end
%Indicare percorso di scrittura file dati x=tempo di integrazione e y=Allan deviation
fileName=sprintf('s%sd%s-%s-%s-%s',serie,disco,power,amplit,Mode);
dlmwrite(percorso+"\"+fileName+".DAT",[tau Allan])

figure(2)
for i=1:n_modi
loglog(tau,Allan(:,i),'o','DisplayName', ['RBM', num2str(i)]);
legend
hold on
end
xlabel('time (sec)');
ylabel('\sigma_{Allan}(\tau)')
set(gca,'fontsize',16)
title(fileName,'fontsize',12)

%indicare percorso scrittura file immagine .png e .fig
box on
print(percorso+"\"+fileName,'-dpng')
savefig(percorso+"\"+fileName+".fig")