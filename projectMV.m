%%%% Ground Truth   %%%% Majority Voting
clc;close all;clear all
Fs = 2000;
T =1/Fs;
N = 1000;
gap = 20;              %%%eleting Posture Between Motions
motion = 1; % motion az 1 ta 9
Tekrar = 1; % Tekrar az 1 ta 10
nstate = 10*(motion-1)+(Tekrar);
w =hamming(N)';
load('C:\Users\Microsoft\Desktop\JodasazaiHarekatVaziat\Data\S1_E1_A1.mat')
DCH = {[2 3 6],[1 2],[5 6],[7 8],[11 12 13 14],[8 11 14],[2 6 8 11 14],[2 5 6 7 8 10 11 13 14],[6 8 11 14]};
%x = emg(:,1);
%[b,a]=butter(2,5/Fs,'low');
%x = sqrt(x.^2);
%x = 1e5*filter(b,a,x);
l = repetition;
l = diff(l);
p = find(l>0);
q = find(l<0);  %%%

m = 50;  %Down Sampling

%x = x(1 : m : end);
l = l(1 : m : end);

Fs = Fs/m;
t1 = [0:length(l)-1]/Fs;
t2 = t1 ; % n/(m*Fs);

%x = x*5;
l = l*5;


p = round(p/m);
n = round(q/m);
flag = 1;
LV = [];
mv=0;
md=0;
K = 1;

ThFlag = 2;           %%########## 0 & 1 & 2


for k = DCH{motion}
    d = glove(:,k);
    b = [-w(1) w(1:end-1)-w(2:end) w(end)]/(T*sum(w));
    dd = filter(b,1,d);
    v = abs(dd)';
    v = v(1 : m : end);
    d = d(1 : m : end);  %%%%%%%%%%%
    d = d+45;
    xv = v([p(nstate):n(nstate)]);
    xd = d([p(nstate):n(nstate)]);
    xv = [xv(fix(N/(2*m)):end) zeros(1,fix(N/(2*m))-1) ];

    %***********************************************

    if ThFlag==0
        Th(k)=mean(abs(xv)); %%%%%%%%%%%
    elseif ThFlag==1
         Th(k)=mean(abs(xv(xv>=2))); 

    else
        x = sort(xv);
        x = x(x<=0.4*max(x));
        Th(k)=mean(x);
    end

    mv=mv+xv;
    md=md+xd;
    lv = 0*xv+1;   %%%%% Selection of Label 1 (Posture) and 2 (Motion)
    lv(xv>=Th(k)) = 2; 
    LV(K,:) = lv;
    K = K+1;
end
for ia = 1:size(LV,2)
    lv = LV(:,ia);
    for ib = 1:2
        nu(ib) = length(find(lv==ib));
    end
    [Pnu,Inu(ia)] = max(nu);
end
for r = [2 1]     %%Deleting the Posture Part Between the Motion Parts 
    y = find(Inu==r);
    I = find(diff(y)>1)
    Is = y(I)               
    Ie = y(I+1)            
    des = Ie-Is;             
    for k = 1:length(Is)
        if des(k)<gap
            Inu(Is(k):Ie(k)) = r;
        end
    end
end
for k = DCH{motion}
    d = glove(:,k);
    b = [-w(1) w(1:end-1)-w(2:end) w(end)]/(T*sum(w));
    dd = filter(b,1,d);
    v = abs(dd)';
    d = d(1 : m : end);
    v = v(1 : m : end);
    d = d+45;
    subplot(4,3,flag)
    %xx = x([p(nstate):n(nstate)]);
    xt2 = t2([p(nstate):n(nstate)]);
    xt1 = t1([p(nstate):n(nstate)]);
    xd = d([p(nstate):n(nstate)]);
    xv = v([p(nstate):n(nstate)]);
    xl = l([p(nstate):n(nstate)]);
    xv = [xv(fix(N/(2*m)):end) zeros(1,fix(N/(2*m))-1) ];
    plot (xt1,xd, 'm');hold on
    plot (xt2,xv, 'y');hold on
    for i=1:length(xd)-1        %%%%%%%%%%%%
        if Inu(i)==2    %%%% Based on Majority Voting
             %if xv(i)>=Th(k)  %%% Based on Threshold
            plot(xt1([i i+1]),xd([i i+1]),'g');
            plot(xt2([i i+1]),xv([i i+1]),'r');
        end
    end
    flag=flag+1;
    plot([xt1(1) xt1(end)],[Th(k) Th(k)],':k')
    xlabel('Time (s)')
    ylabel('Amplitude')
    xlim([xt1(1) xt1(end)])
    title(['Channel ',num2str(k)])
end
subplot(4,3,flag)

%%%length(DCH{motion})=Number of Desired Channels
mv=mv/length(DCH{motion});

hold on
plot (xt2,mv, 'y')

if ThFlag==0
    Th=mean(abs(mv)); %%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&
elseif ThFlag==1
         Th=mean(abs(mv(mv>=2)));
else
        x = sort(mv);
        x = x(x<=0.4*max(x));
        Th=mean(x);
end

for i=1:length(xd)-1
    % if Inu(i)==2           %%%% Based on Majority Voting
    if mv(i)>=Th         
        plot(xt2([i i+1]),mv([i i+1]),'r');
    end
end
plot([xt1(1) xt1(end)],[Th Th],':k')
xlabel('Time (s)')
ylabel('Amplitude')
xlim([xt1(1) xt1(end)])
title('Mean ')

figure
subplot(4,1,1)
plot(xt2,LV(1,:),'r')
xlabel('Time (s)')
ylabel('Label')
title('Channel 2 ')

subplot(4,1,2)
plot(xt2,LV(2,:),'g')
xlabel('Time (s)')
ylabel('Label')
title('Channel 3 ')

subplot(4,1,3)
plot(xt2,LV(3,:),'m')
xlabel('Time (s)')
ylabel('Label')
title('Channel 6 ')

subplot(4,1,4)
plot(xt2,Inu,'k')
xlabel('Time (s)')
ylabel('Label')
title('Majority Voting')
