function Uscof=ufem_vs(uvp_on,fs,ks,wo,count)
[nF,a1,a1g,a1g0,a10]=deal(uvp_on{1:5,1});
[Uscof]=ups_online1(wo',nF,a1,a1g,a1g0,a10,fs,ks);
for i=2:count
    [nF,a1,a1g,a1g0,a10,Aq,Aq0]=deal(uvp_on{:,i});
    [uvtem]=ups_online(wo',Uscof,nF,a1,a1g,a1g0,a10,Aq,Aq0,fs,ks);
    Uscof=[Uscof uvtem];
end
