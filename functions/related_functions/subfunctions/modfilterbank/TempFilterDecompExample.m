[X,Fs]=audioread('/Users/escabi/data/DanSanes/SanesLabSounds/TextureMaskingParadigm/foredirectoryspeech/Speech_Trial001.wav')
[CochData]=cochleogram(X,Fs,1/10,100,10000,500,2,'Amp');
 
S=CochData.SdB;
Fst=1/CochData.taxis(2);
X=log2(CochData.faxis/CochData.faxis(1));
Fss=1/X(2);
 
 
[ModBankParam] = ModFiltBankParamGaborAlpha1ns([1 1 16 512 .2 3.2 .5 .5]);
STRFBank = STRFModFiltBankGaborAlpha1ns(ModBankParam,[Fst Fss]);
Lt=size(S,2);
Ls=size(S,1);
 
MaxNt=-9999;
for k=1:size(STRFBank.F,1)
    for l=1:size(STRFBank.F,2)
   
            MaxNt=max(size(STRFBank.F(k,l).H,2),MaxNt);
            
    end
end
YY=zeros(Ls,Lt+MaxNt-1,size(STRFBank.F,1),size(STRFBank.F,2));
 
for k=1:size(STRFBank.F,1)
    for l=1:size(STRFBank.F,2)
 
            
            Y=conv2(S,STRFBank.F(k,l).H);
            Lt=size(S,2);
            Ls=size(S,1);
            
            Ns=(size(STRFBank.F(k,l).H,1)-1)/2;
            Nt=size(STRFBank.F(k,l).H,2);
            
            YY(:,1:Lt+Nt-1,k,l)=Y(Ns+1:end-Ns,:);
 
            
            %SS=SS(:,1:Lt+Nt-1);
%             
%             size(S)
%             size(SS)
%             size(STRFBank.F(k,l).H)
%             
%             
%             imagesc(SS)
%             pause
            
        
        [k l]
    end
end
 
 
Max=max(max(max(max(YY))));
count=1
for k=1:2:18
    for l=1:2:43
        
        
        subplot(9,22,count)
        imagesc(YY(:,:,k,l))
        caxis([-Max Max]*.1)
        count=count+1
    end
end
 
for k=1:18
    for l=1:43
        MP(l,k)=sum(sum(YY(:,:,k,l).^2));
    end
end 

for k=1:18
    for l=1:43
        MP(l,k)=var(reshape(YY(:,:,k,l),1,numel(YY(:,:,k,l))));
    end
end
 
 