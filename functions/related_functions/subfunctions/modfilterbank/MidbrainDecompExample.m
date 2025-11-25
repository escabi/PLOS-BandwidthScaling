[X,Fs]=audioread('/Users/escabi/data/DanSanes/SanesLabSounds/TextureMaskingParadigm/foredirectoryspeech/Speech_Trial001.wav');
[CochData]=cochleogram(X,Fs,1/10,100,10000,500,2,'Amp','n',60,'GammaTone','erb','hil','bspline');
 
beta=[1 1 4 512 .1 3.2 1 1];
[MidData]=midbrainogram(CochData,beta)

N1=size(MidData.Y,3);
N2=size(MidData.Y,4);
Max=max(max(max(max(MidData.Y))));
for k=1:N1
    for l=1:N2

        subplot(N1,N2,l+(k-1)*N2)
        imagesc(MidData.Y(:,:,k,l)), colormap jet, set(gca,'YDir','normal')
        %set(gca,'visible','off')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        %caxis([-Max Max])
        if k==1
            title(['Fm=' num2str(MidData.ModBankParam.F(k,l).Fm,3) ])        
        end
        if l==1
            ylabel(['RD=' num2str(MidData.ModBankParam.F(k,l).RD,3) ]) 
        end
        
    end
end


%Putting through Sigle Layer Spiking Network
L=size(MidData.Y,1);
Tau=5;
Tref=1;
Nsig=3;
SNR=-5;
SigE=0.0255;
SigI=1.5*SigE;
EIR=1.5
Fs=1/CochData.taxis(2);
flag=3;
detrendim='n'
detrendin='n'

NF=size(CochData.SdB,1);
NT=size(CochData.SdB,2);
for k=1:N1
    for l=1:N2
        X=MidData.Y(:,:,k,l);
        [Y]=integratefirenetworkcontmulti(X,L,Tau,Tref,Nsig,SNR,SigE,SigI,EIR,Fs,flag,detrendim,detrendin);
        MidData.S(:,:,k,l)=Y;
    end
end


for k=1:N1
    for l=1:N2
         subplot(N1,N2,l+(k-1)*N2)
        [i,j]=find(MidData.S(:,:,k,l));
        plot(j/Fs,i,'k.')
        ylim([0 NF])
        xlim([0 NT/Fs])
        %set(gca,'visible','off')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
 
        if k==1
            title(['Fm=' num2str(MidData.ModBankParam.F(k,l).Fm,3) ])        
        end
        if l==1
            ylabel(['RD=' num2str(MidData.ModBankParam.F(k,l).RD,3) ]) 
        end
        
    end
end


% for k=1:18
%     for l=1:43
%         MP(l,k)=sum(sum(YY(:,:,k,l).^2));
%     end
% end 
% 
% 
% for k=1:size(YY,3)
%     for l=1:size(YY,4)
%         MP(l,k)=var(reshape(YY(:,:,k,l),1,numel(YY(:,:,k,l))));
%     end
% end