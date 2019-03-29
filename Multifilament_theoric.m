%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Plasma current centroid position reconstruction%%%%%%%%
%%%%%% Multifilaments,7 filaments , SVD matrix %%%%%%%%%%
tic
close all
clear all
%%% Load Shot
load('shot_46104.mat');
time=1e-3*data.time; %%%% time in ms


%%% Draw the vessel
th = 0:pi/50:2*pi;
xvess = 9 * cos(th)+46;
yvess = 9 * sin(th) ;

%%% Mirnov positions
ang=-15;
for i=1:12
R_mirn(i)=9.35*cosd(ang)+46;
z_mirn(i)=9.35*sind(ang);
ang=ang-30;
end

%%%%%% Lets draw the plasma filaments
th1 = 0:pi/50:2*pi;
R_pls=46;
z_plsm= 0;


degr=30;
radius=5.5; %%% in [cm] (distance from the center of the chamber to the filaments)
nfil=12; %%% Number of filaments
deg_fact=360/(nfil);

for i=1:nfil
    R_filaments(i)=(46)+radius*cosd(degr);
    z_filaments(i)=radius*sind(degr);
    degr=degr+deg_fact;
end



%%%Experimental mesurements[Wb]

%Mirnv_10_fact=1.2803;
time_ins=160;
time_index=find(time == time_ins); %%% Select a time moment where there is plasma current! in [ms]

%%%%%%%%%% Find the exprimental values for that time moment

%%%%without external flux correction

Mirnv_flux(:)=data.mirnv_corr(:,time_index);
%Mirnv_flux(10)=Mirnv_10_fact*Mirnv_flux(10);

Mirnv_flux_corr(:)=data.mirnv_corr_flux(:,time_index);
%Mirnv_flux_corr(10)=Mirnv_10_fact*Mirnv_flux_corr(10);

%%%%% Let's go from [Wb] to {T]
%%%% fmincon needs data to be double
Mirnv_B_exp=double(Mirnv_flux/(50*49e-6)); %%%% [T]
Mirnv_B_exp_corr=double(Mirnv_flux_corr/(50*49e-6)); %%%% [T]

%%%%% Minimization 7 degrees of freedom
%   fval_multi_corr=fminsearch(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,z_filaments(1),R_filaments(1),x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),R_filaments,z_filaments,R_mirn,z_mirn,nfil),[500,500,500,500,500,500,500,500,500,500,500])
%   fval_multi_corr=fminsearch(@(x) ErrorMirnFuncMultiFilam(Mirnv_B_exp_corr,z_filaments(1),R_filaments(1),x,R_filaments,z_filaments,R_mirn,z_mirn,nfil),[500,500,500,500,500,500,500,500,500,500,500])


 %%
% xx_multi_corr=BmagnMultiModule_correct(z_filaments(1),R_filaments(1),fval_multi_corr,R_filaments,z_filaments,R_mirn,z_mirn,nfil);
%%%% Matrix whose elements gives the contribution  to the measuremnt i  to
%%%% a unitary current in the filament j [T]
for i=1:12
    for j=1:nfil
   
         Mfp(i,j)=Bmagnmirnv(z_filaments(j),R_filaments(j),1,R_mirn(i),z_mirn(i)) ;
    end
end

Mpf=pinv(Mfp);
% Mpf=inv(Mfp);
I_filament=Mpf*(Mirnv_B_exp_corr');
I_filament_uncrr=Mpf*(Mirnv_B_exp');

%% Calculate Biot-Savart with the current values from SVD decomposition
xx_multi_SVD=BmagnMultiModule_correct(I_filament,R_filaments,z_filaments,R_mirn,z_mirn,nfil);
xx_multi_SVD_uncrr=BmagnMultiModule_correct(I_filament_uncrr,R_filaments,z_filaments,R_mirn,z_mirn,nfil);

% for i=1:12
%    xx_multi_theo(i) =0;
%     for j=1:7
%         
% xx_multi_theo(i)=Bmagnmirnv(z_filaments(j),R_filaments(j),I_filament(j),R_mirn(i),z_mirn(i)) +xx_multi_theo(i);
%     end
% end
%% Error
RMSE_optim_theo=sqrt(mean((xx_multi_SVD(:)-Mirnv_B_exp_corr(:)).^2));



%% Compute Centroid position
I_filament_all=Mpf*(data.mirnv_corr_flux)/(49*50*1e-6);
I_filament_all_uncrr=Mpf*(data.mirnv_corr)/(49*50*1e-6);
for(i=1:length(I_filament_all))
z0(i)=0.01*sum((z_filaments'.*I_filament_all(:,i)))./sum(I_filament_all(:,i));
r0(i)=0.01*sqrt(sum((R_filaments'.^2).*I_filament_all(:,i))./sum(I_filament_all(:,i)))-0.46;
sumIfil(i)=sum(I_filament_all(:,i));
sumIfil_uncrr(i)=sum(I_filament_all_uncrr(:,i));
end

for(i=1:length(z0))
if(imag(z0(i)~=0))
z0(i)=0.085;
end
if(z0(i)>0.085||z0(i) <-0.085)
z0(i)=0.085;
end
end

for(i=1:length(r0))
if(imag(r0(i))~=0)
r0(i)=0.085;
end
if(r0(i)>0.085||r0(i) <-0.085)
r0(i)=0.085;
end
end


%% Plotting
close all
figure(9)
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mirnv_B_exp_corr ,'-o')
hold on
plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_SVD,'-s')
% plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_SVD_uncrr,'-*')
% plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*xx_multi_corr,'-*')
% plot([1,2,3,4,5,6,7,8,9,10,11,12],1000*Mfp*I_filament,'-s')
grid on
title(['Shot #46104  t= ',num2str(time_ins), '[ms]  Magnetic Probes'])
legend('Experimental Data ',' SVD reconstruction')
xlabel('Mirnov #')
ylabel('Optimization [mT]')
axis equal
%%
close all
figure(11)
plot(time,data.mirnv_corr_flux(1,:))
hold on
plot(time,data.mirnv_SVD_recons(1,:))
xlabel('Time [ms]')
ylabel('Magnetic flux [Wb]')
title(['Shot #46104   Magnetic Probe #1'])
legend('Mirnv #1 flux corrected','Mirnv #1 SVD reconstruction ')
xlim([0 1100])
grid on
%%
close all
figure(16)
index1=find(time == 205.5)
index2=find(time == 258.1)
suptitle('Shot # 46104')
subplot(3,1,1)
plot(time(index1:index2),data.Ip_magn(index1:index2))
ylabel('Current [A]')
title('Plasma current')
grid on
subplot(3,1,2)
plot(time(index1:index2),data.R0(index1:index2))
ylabel('Position [m]')
title('Radial Centroid Position')
grid on
subplot(3,1,3)
plot(time(index1:index2),data.z0(index1:index2))
xlabel('Time [ms]')
ylabel('Position [m]')
title('Vertical Centroid Position')
grid on
%%
%%
close all
figure(10)
plot(time,sumIfil)
hold on
% plot(time,sumIfil_uncrr)
plot(time,data.Ip_magn)
grid on
xlabel('Time [ms]')
ylabel('Current [A]')
legend('\Sigma(I_{fil})','Plasma current')
title('Shot #45680 Plasma Current')

%%
figure(3)
plot(xvess,yvess,'k','linewidth',2)
hold on
plot(46,0,'.m','MarkerSize',790)
plot(R_mirn,z_mirn,'sk','MarkerSize',17)

% plot(fval_multi_corr(2),fval_multi_corr(1),'.k','MarkerSize',20)
for i=1:nfil
    plot(R_filaments(i),z_filaments(i),'.b','MarkerSize',20)
end
    for i = 1:12
    text(R_mirn(i),z_mirn(i),num2str(i),'Color','r','FontSize',13) 


end
% plot([35,57],[0,0])
% plot([46,46],[-9,9],'k')
text(57,0,'LFS','FontSize',15)
text(33,0,'HFS','FontSize',15)
ylim([-11,11])
xlabel('R[cm]')
ylabel('Z[cm]')
grid on
axis equal
toc
<<<<<<< HEAD
%%

=======

%%
>>>>>>> a1512301c4942a96f47983a7e144a2c85859ef28

figure(1)
plot(time,r0)
grid on
xlabel('Time [ms]')
ylabel('Position [m]')
title('Radial centroid position')

figure(2)
plot(time,z0)
grid on
xlabel('Time [ms]')
ylabel('Position [m]')
title('Vertical centroid position')