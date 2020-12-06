close all
clear all


%% Read signal


snrfile_final=zeros(10,9);
i_file=11;
k_snr=1;
while(i_file<=20)

filename = (['a00',num2str(i_file),'.wav']);

[x_original, fs] = audioread(filename);
x_original=x_original(:,1);
m=max(abs(x_original));
x_original = x_original/max(abs(x_original)); % normalize signal


%% Clip signal:  %% first no clip libin

ssinput=[3.5,7,10.8,14,18,22,27,33,42];% speech.wav

ii=9;
snr_final=zeros(1,9);

while(ii>=1)
SNRInput = ssinput(ii);


[y_clip, ~] = clip_signal(x_original, SNRInput);

z=y_clip; %%原始信号，剪裁长度
max_temp=max(abs(z));

SNRin = SNR(x_original,y_clip);
fprintf('Input SNR: %.3f dB\n',SNRin)

% % Detect clipping level:
ClippingLevel = max(abs(y_clip));
clipped = (abs(y_clip) - ClippingLevel).^2 < 1e-4;
reliable_samples = ~clipped;

SNRin_clipped = SNR(x_original(~reliable_samples),y_clip(~reliable_samples));
fprintf('%.1f percent of clipped samples\n', sum(~reliable_samples)/length(x_original)*100)

frame_len=1024;

xori_len=length(x_original);
k=ceil(xori_len/frame_len);

x_add=zeros(k*frame_len,1);
x_add(1:xori_len,1)=x_original;
x_matrix=reshape(x_add,[frame_len,k]);

y_add=zeros(k*frame_len,1);
y_add(1:xori_len,1)=y_clip;
y_matrix=reshape(y_add,[frame_len,k]);

z_add=zeros(k*frame_len,1);
z_matrix=reshape(z_add,[frame_len,k]);


count_frame=1;
while (count_frame<=k)

x=x_matrix(1:frame_len,count_frame);
y=y_matrix(1:frame_len,count_frame);

Nsize=frame_len;%4096,16384,65536,length(x)


%% Clip signal: 

b=1;
j=1136;
while(j == 136)
i=1;
y_max=max(abs(y));    
if y_max < max_temp
    i=2000;
end;
    
d = 64; % signal dimension 64
l = j; % co-sparsity       40
X = im2col(y, [d, 1], 'sliding');% attention 20200219
ZZ=X;

reliable_samples_mat= (X<ClippingLevel & X>-ClippingLevel);
clipped_pos_mat = (X==ClippingLevel);
clipped_neg_mat = (X==-ClippingLevel);

Xntrain=X;

p=144;

% OmegaInit = normr(randn(p, d));

KK=144;
Pn=ceil(sqrt(KK));  
bb=8;
DCT=zeros(bb,Pn); for kk=0:1:Pn-1,         
    V=cos([0:1:bb-1]'*kk*pi/Pn);    
    if kk>0, V=V-mean(V); 
    end;    
    DCT(:,kk+1)=V/norm(V);
end; 
G=kron(DCT,DCT);
OmegaInit=G';
% OmegaInit = normr(randn(p, d));



paramASimCO.initialDictionary = OmegaInit;
temp_Omega=OmegaInit;

paramASimCO.itN =1;
paramASimCO.cosparsity = l;
paramASimCO.numOmegaUpdate = 1;


lambda = 0.1;

temp3=NaN;


cost_norm1 = NaN(200000,1); % save cost at each iteration
cost_norm2 = NaN(200000,1); % save cost at each iteration.
snrr = NaN(200000,1); % save cost at each iteration.
temp21=NaN;
temp22=100000000;
temp11=NaN;
temp12=NaN;
y_max=max(abs(y));    
if y_max < max_temp
    i=2000;
end;
snrr_temp=0;
while (i <=2000)

Xntrain_rdc = double(Xntrain);

flag11=1;
i1=1;
% decreasing

 while(flag11 == 1) 
      
[Omega,f_itn,X_sparse] = normalizedRowOmega_sad(Xntrain_rdc, paramASimCO);

if ((f_itn-temp11 > 0) | (f_itn-temp12 > 0));
    flag11=1;
else
    flag11=0;    
end;
cost_norm2(i)=norm(Omega-temp_Omega,"fro");
temp_Omega=Omega;
paramASimCO.initialDictionary = Omega;
 temp11=f_itn;
i1=i1+1;
end;
temp12=f_itn;

Y=Xntrain;
Ym = mean(Y);


 flag21=1;
 i2=1;

 while(i2 <=1) 
    
    
[Xdn,cost] = ADxh108_reliable(Y, Omega,reliable_samples_mat,X_sparse,ClippingLevel,clipped_pos_mat,clipped_neg_mat,ZZ);

 i2=i2+1;
end;

Xntrain = Xdn ;
 
if (rem(i,2000)==0) 
 Xdn=Xntrain;
% Reconstruct
[MM,NN]=size(y);
Idn=col2imstep(Xdn,[MM NN],[d 1]);
cnt=countcover([MM NN],[d  1],[1 1]);



Idn=Idn./cnt;
Idn=reshape(Idn,Nsize(1),1);
Idnn=Idn;
counti=round(i/1);
snrr(counti)=SNR(x,Idnn);
fprintf(' SNRin %.2f. SNRout %.4f. File is a00%.0f.\n',SNRInput,snrr(i),i_file)  
end ;
 cost_norm1(1)=1000000000000;
 cost_norm1(i+1)=cost+f_itn;
 cost_temp=abs(cost_norm1(i+1)-cost_norm1(i)); 
  

 i=i+1; 
end
j=j+8;
a(b,:)=snrr;
b=b+1;
end;
z_matrix(1:frame_len,count_frame)=Idnn;
count_frame=count_frame+1;
end;

z_temp=reshape(z_matrix,[],1);
z=z_temp(1:xori_len,1);

snr_final(1,ii)=SNR(x_original,z);
% ii=ii+1;
ii=ii-1;
end;
snrfile_final(k_snr,1:9)=snr_final;
save(['a00',num2str(i_file),'.mat']);
i_file=i_file+1;
k_snr=k_snr+1;
end;

SNRSAD=mean(snrfile_final);

p=0.1:0.1:0.9;
figure,plot(p,SNRSAD,'--ko','LineWidth',1);
set(gca,'linewidth',0.5,'fontsize',15);
axis([0.1 0.9,0 60])
grid on;
xlabel('Clipping level');
ylabel('SNR');
legend('SAD');
