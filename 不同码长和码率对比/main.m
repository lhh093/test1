clear all;
close all;
clc;

%码率R=(n-m)/n
% %码长相同，不同码率（R=1/2，N=1024）
% m1=512;
% n1=1024;
% [someSNR_1, ave_BER50_1] = BER(m1,n1);

%%码长相同，不同码率（R=2/3，N=1024）
% m2=341;
% n2=1024;
% [someSNR_2, ave_BER50_2] = BER(m2,n2);

% %码长相同，不同码率（R=3/4，N=1024）
% m3=256;
% n3=1024;
% [someSNR_3, ave_BER50_3] = BER(m3,n3);


%码率相同，不同码长（N=1024，R=1/2）
m1=512;
n1=1024;
[someSNR_1, ave_BER50_1] = BER(m1,n1);

%码率相同，不同码长（N=4096，R=1/2）
m2=1024;
n2=2048;
[someSNR_2, ave_BER50_2] = BER(m2,n2);

%码率相同，不同码长（N=16374，R=1/2）
m3=2048;
n3=4096;
[someSNR_3, ave_BER50_3] = BER(m3,n3);
figure
plot(someSNR_1,ave_BER50_1,'g-*');
hold on;
plot(someSNR_2,ave_BER50_2,'r-+');
hold on;
plot(someSNR_3,ave_BER50_3,'k-o');
hold on;
set(gca,'Yscale','log');
ylabel('BER');
xlabel('SNR');
grid on;
%title('不同码率对比');
 title('不同码长对比');
%legend('R=1/2','R=2/3','R=3/4');
legend('N=256','N=4096','N=16374');