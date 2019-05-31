figure()

plot(ghVdsOutTTFS(1,:,2),ghVdsOutTTFS(1,:,3),'-k',ghVdsOutTTFS(2,:,2),ghVdsOutTTFS(2,:,3),'.-r', ...
    ghVdsOutTTFS(3,:,2),ghVdsOutTTFS(3,:,3),':k',ghVdsOutTTFS(4,:,2),ghVdsOutTTFS(4,:,3),'-r', ...
    ghVdsOutTTFS(5,:,2),ghVdsOutTTFS(5,:,3),':r','LineWidth',2)

figure()
symbolbag={'-ks','-k^','-ko',':ks',':k^',':ko','-kp'};
for z=1:size(Lippertghvals,2)
    plot(ghVdsOutTTFS(z,:,2),ghVdsOutTTFS(z,:,3),symbolbag{z})
    hold on;
end
h=legend(strlegend)
xlabel('V_{ds}^{out} (mV)')
ylabel('TTFS (ms)')