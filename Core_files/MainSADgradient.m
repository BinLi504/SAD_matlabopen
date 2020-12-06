close all
clear all



%% Read signal


filename = 'dev_male2_150ms_1_ch12.wav';
 
[x, fs] = audioread(filename);

x=x(1:16384,1);
Nsize=length(x);

m=max(abs(x));
x = x/max(abs(x)); % normalize signal


%% Clip signal:

thans=zeros(1,10);
a=zeros(100,200000);
bb=zeros(100,5000);

SNRInput =13.4;
[y_clip, ~] = clip_signal(x, SNRInput);
z=y_clip;

SNRin = SNR(x,y_clip);
fprintf('Input SNR: %.3f dB\n',SNRin)

% % Detect clipping level:
ClippingLevel = max(abs(y_clip));
clipped = (abs(y_clip) - ClippingLevel).^2 < 1e-4;
reliable_samples = ~clipped;
SNRin_clipped = SNR(x(~reliable_samples),y_clip(~reliable_samples));


d = 64; % signal dimension 64
l = 120; % co-sparsity  
X = im2col(y_clip, [d, 1], 'sliding');% attention 20200219
ZZ=X;

reliable_samples_mat= (X<ClippingLevel & X>-ClippingLevel);
clipped_pos_mat = (X==ClippingLevel);
clipped_neg_mat = (X==-ClippingLevel);

Xntrain=X;

p = 2*d; % number of atoms

OmegaInit = normr(randn(p, d));

paramASimCO.initialDictionary = OmegaInit;


paramASimCO.itN =1;
paramASimCO.cosparsity = l;
paramASimCO.numOmegaUpdate = 1;


lambda = 0.1;


temp3=NaN;

i=1;
%save cost value
f_itn_final = NaN(5000,1); % save cost at each iteration 
cost_final = NaN(5000,1); % save cost at each iteration
cost_norm1 = NaN(5000,1); % save cost at each iteration
cost_norm2 = NaN(5000,1); % save cost at each iteration.
snrr = 0; % save cost at each iteration.
temp21=NaN;
temp22=100000000;
temp11=NaN;
temp12=NaN;
while (i <=5000)
Xntrain_rdc = double(Xntrain);
flag11=1;
i1=1;

while(flag11 == 1) 
  
[Omega,f_itn,X_sparse] = normalizedRowOmega_sad(Xntrain_rdc, paramASimCO);

if ((f_itn-temp11 > 0) | (f_itn-temp12 > 0));
    flag11=1;
else
    flag11=0;    
end;

paramASimCO.initialDictionary = Omega;
temp11=f_itn;
i1=i1+1;
end;
temp12=f_itn;
f_itn_final(i)=f_itn;
Y=Xntrain;

Ym = mean(Y);

 flag21=1;
 i2=1;
  
 while(i2 <=1) 
[Xdn,cost] = SAD_gradient(Y, Omega,reliable_samples_mat,X_sparse,ClippingLevel,clipped_pos_mat,clipped_neg_mat,ZZ);
i2=i2+1;
end;
temp22=cost;   

cost_final(i)=cost;

Xntrain = Xdn ;
 
if (rem(i,5000)==0) 
 Xdn=Xntrain;
% Reconstruct
[MM,NN]=size(y_clip);
Idn=col2imstep(Xdn,[MM NN],[d 1]);
cnt=countcover([MM NN],[d  1],[1 1]);
Idn=Idn./cnt;
Idn=reshape(Idn,Nsize(1),1);
Idnn=Idn;
counti=round(i/1);
snrr=SNR(x,Idnn);
end ;

i=i+1;

end
fprintf('SNRin %.2f. SNRout %.4f. \n',SNRInput,snrr)  
figure, plot(1:Nsize, x, 1:Nsize, Idnn,1:Nsize,z,'--');
legend('clean A','estimate X','clipped Y')

