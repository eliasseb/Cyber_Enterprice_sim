

theta=0:10:360;
ETA=[];
EETA=[];
for i=theta

    Eta=Eta_d_elipse(i)';
    EETA=[EETA Eta];
    
    quiver(Eta(2),Eta(1),sin(Eta(3)),cos(Eta(3)));
    ETA=[ETA Eta(3)];
    hold on

end

figure(2)
plot(ETA)