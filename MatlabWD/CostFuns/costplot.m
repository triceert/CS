function costplot

close all

recapex=[ 0.674 ; 0.736 ; 8.577 ; 20.325];
reopex=[ 5.175 ; 5.201 ; 5.172 ; 14.831];
nh3abscapex=[ 0.661 ; 0.663 ; 0.691 ; 0.816];
nh3absopex=[ 0.136 ; 0.137 ; 0.372 ; 0.146];
hcnabscapex=[ 0.230 ; 0.230 ; 0.236 ; 0.275];
hcnabsopex=[0.040 ; 0.040 ; 0.045 ; 	0.045];
hcndistcapex=[1.272 ; 1.272 ; 1.511 ; 1.536];
hcndistopex=[0.526 ; 0.525 ; 0.663 ; 0.683];

totex=horzcat(recapex+reopex,nh3abscapex+nh3absopex,hcnabscapex+hcnabsopex,hcndistcapex+hcndistopex);
capex=horzcat(recapex,nh3abscapex,hcnabscapex,hcndistcapex);
opex=horzcat(reopex,nh3absopex,hcnabsopex,hcndistopex);

c = categorical({'Case 1','Case 2','Case 3','Case 4'});


figure
hold on
subplot(2,2,3)
bar(c,totex,'stacked')
box on 
pbaspect([1 1 1]);


ylabel('$Totex/[Mio. \$]$')


subplot(2,2,1)
hold on
bar(c,capex,'stacked')
legend({'Reactor','$NH_3$ Absorber','$HCN$ Absorber','Distillation'},'Location','Northwest')
box on 
pbaspect([1 1 1]);

ylabel('$Capex/[Mio. \$]$')

subplot(2,2,2)
hold on
bar(c,opex,'stacked')
box on 
pbaspect([1 1 1]);

ylabel('Opex/[Mio. \$/y]')

end