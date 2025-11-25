function MtfBank = MtfModFiltBankGaborAlpha1ns(ModBankParam,input)
%Generate Modualtion filters
[Ns,Nt] = size(ModBankParam.F);

for i = 1:Nt
    for j = 1:Ns
        Beta = ModBankParam.F(j,i).Beta;
        H = mtfgaboralpha1modelns(Beta,input);

        %Save Filter
        MtfBank.F(j,i).Beta = Beta;
        MtfBank.F(j,i).input = input;
        MtfBank.F(j,i).H = H;
        clear H;
    end
end

end
