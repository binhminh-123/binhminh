clc;
clear all
close all

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1.Fun. for the MinRate&SumRate_versu_K in near-field.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%最小用户速率随用户数量增长的变化
%% Thiết lập thông số mô phỏng
N1 =64;
N2 = 8;
N = N1*N2; % số phần tử mảng RIS
d = 0.5; % khoảng cách giữa các phần tử lambda/2

%K số người dùng
% K=4;
% U=K;
K_max=20;
K_list=[2:2:K_max];

% ITER là số lần lặp để tính trung bình, trong bài báo này thì ITER = 600;
ITER = 30; lấy ít cho ra số lẹ lẹ

A = 4;
delta = 0.25;
D_oversample=1;

realsnr=5;
SNR_linear = 10.^(realsnr/10.);
%% Thiết lập codebook
disp("Gene Near and Far Codebooks……")
% Thiết lập codebook viễn trường
% Tính theo công thức 5a và 5b
UN1 = exp(1i*2*pi*[0:(N1-1)*D_oversample]'*d*[0:N1-1]*(2/N1/D_oversample))/sqrt(N1);
UN2 = exp(1i*2*pi*[0:(N2-1)*D_oversample]'*d*[0:N2-1]*(2/N2/D_oversample))/sqrt(N2);
far_codebook = kron(UN1,UN2);

% Tọa độ phân bố của BS và user
P3 = [2500*d,-2500*d,1200*d,200*d,0*d,-1000*d];
P4=P3;
P1=P3;
P2=P3;
Delta = 1*[100*d,100*d,100*d,100*d,100*d,100*d];
Delta1 = Delta*A;
[near_codebook1,record] = generate_near_field_codebook(N1,N2,d,P3,P4,Delta1);
near_codebook1=near_codebook1./sqrt(N);
disp("Finish Codebooks Gene")

%% Lưu số liệu tốc độ truyền ( để vẽ đồ thị )
Record_SumR_FF_RIS=zeros(length(K_list),1);
Record_SumR_NF_RIS=zeros(length(K_list),1);
Record_SumR_FF_AP=zeros(length(K_list),1);
Record_SumR_NF_AP=zeros(length(K_list),1);
Record_SumR_MM=zeros(length(K_list),1);
Record_SumR_Rand=zeros(length(K_list),1);

Record_MinR_FF_RIS=zeros(length(K_list),1);
Record_MinR_NF_RIS=zeros(length(K_list),1);
Record_MinR_FF_AP=zeros(length(K_list),1);
Record_MinR_NF_AP=zeros(length(K_list),1);
Record_MinR_MM=zeros(length(K_list),1);
Record_MinR_Rand=zeros(length(K_list),1);

%% Hàm chính
t0 = clock;
for idx_K=1:length(K_list)
    LengthK_list=length(K_list);
    num_K=K_list(idx_K);
    K=num_K;
    %% Lưu biến tạm để tính trung bình (a=a+data./ITER)
    temp_SumR_FF_RIS=0;
    temp_SumR_NF_RIS=0;
    temp_SumR_FF_AP=0;
    temp_SumR_NF_AP=0;
    temp_SumR_MM=0;
    temp_SumR_Rand=0;
    
    temp_MinR_FF_RIS=0;
    temp_MinR_NF_RIS=0;
    temp_MinR_FF_AP=0;
    temp_MinR_NF_AP=0;
    temp_MinR_MM=0;
    temp_MinR_Rand=0;
    %% Vòng lặp chính
    for idx_iter=1:ITER
        
        %% Tạo kênh BS tới RIS
        FCCodewordsBuffer=zeros(N,num_K);
        NCCodewordsBuffer=zeros(N,num_K);
        PftCodewordsBuffer=zeros(N,num_K);
        FCGainBuffer=zeros(num_K,1);
        NCGainBuffer=zeros(num_K,1);
        PftGainBuffer=zeros(num_K,1);
        max_index_2ndlayer=zeros(num_K,1);
        % Tạo kênh truyền BS tới RIS
        [G,px1,py1,pz1,alpha] = generate_G_near_field_channel(N1,N2,P1);
        GG=zeros(N,num_K);
        %% Lưu thông số thời gian
        fprintf('For UserNum  (NearField):i_num=%d of %d,i_iter=%d of %d | run %.4f s\n',idx_K,LengthK_list,idx_iter,ITER,  etime(clock, t0));
        
        %% Tạo kênh truyền RIS tới người dùng
        for k=1:num_K
            [hK,px2,py2,pz2,alpha] = generate_hr_near_field_channel(N1,N2,1,P2);
            %hK=hK./sqrt(N);
            Hc = diag(hK)*G;
            GG(:,k)=Hc;           
           %% Tạo tia cận trường và viễn trường：(Get N(or F)CCodewordsBuffer、N(or F)CGainBuffer、PftGainBuffer)
            % Tạo tia viễn trường
            [maxGainFC,idxFC]=max(abs(far_codebook*Hc).^2);
            FCCodewordsBuffer(:,k)=far_codebook(idxFC,:).';
            FCGainBuffer(k)=maxGainFC;
            
            % Tạo tia cận trường
            array_gain = 0;
            max_index=-1;
            for i =1:size(near_codebook1,1)
                if abs(near_codebook1(i,:)*Hc)^2>array_gain
                    max_index=i;
                    array_gain=abs(near_codebook1(i,:)*Hc)^2;
                end
            end
            NCCodewordsBuffer(:,k)=near_codebook1(max_index,:).';
            % Tạo mã 2 lớp
            P21=[record(max_index,1)+Delta1(1)/2,record(max_index,1)-Delta1(1)/2,record(max_index,2)+Delta1(2)/2,record(max_index,2)-Delta1(2)/2,record(max_index,3)+Delta1(3)/2,record(max_index,3)-Delta1(3)/2];
            P22=[record(max_index,4)+Delta1(4)/2,record(max_index,4)-Delta1(4)/2,record(max_index,5)+Delta1(5)/2,record(max_index,5)-Delta1(5)/2,record(max_index,6)+Delta1(6)/2,record(max_index,6)-Delta1(6)/2];
            
            near_codebook2 = generate_near_field_codebook(N1,N2,d,P21,P22,Delta1*delta);
            near_codebook2=near_codebook2./sqrt(N);
            
            for i =1:size(near_codebook2,1)
                if abs(near_codebook2(i,:)*Hc)^2>array_gain
                    array_gain=abs(near_codebook2(i,:)*Hc)^2;
                    max_index_2ndlayer(k)=i;
                end
            end
            if max_index_2ndlayer(k)>0 %二层码本起了作用
                NCCodewordsBuffer(:,k)=near_codebook2(max_index_2ndlayer(k),:).';
                maxGainNC=array_gain;
                NCGainBuffer(k)=abs(near_codebook2(max_index_2ndlayer(k),:)*Hc)^2;
            else
                
                NCGainBuffer(k)=abs(NCCodewordsBuffer(:,k).'*Hc)^2;
                
            end%
            
            %Pft precoding
            wc_opt = exp(1j*phase(Hc'));%If can not work in this line maybe with different Matlab version, you can try to use "angle" function to replace "phase"
            wc_opt=wc_opt./abs(wc_opt)/sqrt(N);
            array_gainpft = abs(wc_opt*Hc)^2;
            PftCodewordsBuffer(:,k)=wc_opt.';
            PftGainBuffer(k)=array_gainpft;
        end
        %% Xử lý độ lợi kênh truyền
        Product_mxg_DFT=prod(sqrt(FCGainBuffer));
        MultiBeamFC_Orig=FCCodewordsBuffer*((Product_mxg_DFT./sqrt(FCGainBuffer)));
        Product_mxg_NC=prod(sqrt(NCGainBuffer));
        MultiBeamNC_Orig=NCCodewordsBuffer*((Product_mxg_NC./sqrt(NCGainBuffer)));
        Product_mxg_Pft=prod(sqrt(PftGainBuffer));
        MultiBeamPft_Orig=PftCodewordsBuffer*((Product_mxg_Pft./sqrt(PftGainBuffer)));
        %% Chồng chập tia viễn trường
        %MultiBeamFC_Orig=sum(FCCodewordsBuffer,2);
        record_zeroFC=find(MultiBeamFC_Orig==0);
        MultiBeamFC_Orig(record_zeroFC)=exp(1j*2*pi*rand)/sqrt(N);
        MultiBeamFCRIS=MultiBeamFC_Orig./abs(MultiBeamFC_Orig)/sqrt(N);
        MultiBeamFCAP=MultiBeamFC_Orig./max(abs(MultiBeamFC_Orig))/sqrt(N);
        %% Chồng chập tia cận trường
        %MultiBeamNC_Orig=sum(NCCodewordsBuffer,2);
        record_zeroNC=find(MultiBeamNC_Orig==0);
        MultiBeamNC_Orig(record_zeroNC)=exp(1j*2*pi*rand)/sqrt(N);
        MultiBeamNCRIS=MultiBeamNC_Orig./abs(MultiBeamNC_Orig)/sqrt(N);
        MultiBeamNCDig=MultiBeamNC_Orig/norm(MultiBeamNC_Orig);%Use for MM's v_abs setting
        MultiBeamNCAP=MultiBeamNC_Orig./max(abs(MultiBeamNC_Orig))/sqrt(N);%Which means that the Amplitude & Phase can be adjusted. (But the Amplitude is in [0,1])
        %
        %% Chồng chập tia Pft ( kết hợp Tiền xử lý và Kết hợp tín hiệu vào tạo tia )
        %MultiBeamPft_Orig=sum(PftCodewordsBuffer,2);
        record_zeroPft=find(MultiBeamPft_Orig==0);
        MultiBeamPft_Orig(record_zeroPft)=exp(1j*2*pi*rand)/sqrt(N);
        MultiBeamPftRIS=MultiBeamPft_Orig./abs(MultiBeamPft_Orig)/sqrt(N);
        MultiBeamPftDig=MultiBeamPft_Orig/norm(MultiBeamPft_Orig);%Use for MM's v_abs setting
        MultiBeamRand=exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N);
        %
        %% Tính độ lợi đa tia
        G_MultiBeam_FFSuperpose=abs(MultiBeamFCRIS.'*GG).^2;%size:1 * K (RIS BF)
        G_MultiBeam_NFSuperpose=abs(MultiBeamNCRIS.'*GG).^2;%(RIS BF)
        G_MultiBeam_PftSuperpose=abs(MultiBeamPftDig.'*GG).^2;%(Dig BF)
        G_MultiBeam_FFSuperpose_AP=abs(MultiBeamFCAP.'*GG).^2;%size:1 * K (AP BF)
        G_MultiBeam_NFSuperpose_AP=abs(MultiBeamNCAP.'*GG).^2;%(AP BF)
        G_MultiBeam_Rand=abs(MultiBeamRand.'*GG).^2;%(AP BF)

        %% Thuật toán MM
        loss3 = [];
        eta=1;
        v_abs3=max(sqrt(G_MultiBeam_PftSuperpose))*ones(1,num_K)'*eta;
        w3 = exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N);
        
        v_phase=exp(1j*1*pi*(2*rand(num_K,1)-1));
        v3 = v_abs3.*v_phase;
        
        A=conj(NCCodewordsBuffer).';%根据近场码本确定A     k*n
        A=diag(sqrt(NCGainBuffer))*conj(NCCodewordsBuffer).';%根据近场码本确定A     k*n
                
        %A=GG.';
        lambda = max(eig(A'*A));
        max_iter=200;
        % Vòng lặp
        for i=1:max_iter
            % update v_phase
            f = A*w3;
            v_phase = exp(1j * angle(f)); % cập nhật biến phụ
            v3 = v_abs3 .* v_phase;
            % cập nhật thông số w3
            w3_0iter=w3;
            for j_MM=1:10
                temp = A' * v3 - A' * A * w3_0iter + lambda * w3_0iter ;
                w3_0iter = exp( 1j*angle( temp ))/sqrt(N);
                
            end
            w3=w3_0iter;
            loss3 = [loss3, norm(v3 - A * w3, 2)];
        end
        %% Tính độ lợi đa tia bằng MM
        G_MultiBeam_MM=abs(w3.'*GG).^2;
        G_MM=abs(w3.'*GG).^2;
        %% Tính tốc độ truyền cho từng người dùng
        R_FF_RIS=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose.');%Size K*1
        R_NF_RIS=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose.');
        R_FF_AP=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose_AP.');
        R_NF_AP=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose_AP.');
        R_MM=log2(1+SNR_linear.*G_MultiBeam_MM.');
        R_Rand=log2(1+SNR_linear.*G_MultiBeam_Rand.');
        %% Tổng tốc độ truyền & Tốc độ truyền tối thiểu
        SumR_FF_RIS=sum(R_FF_RIS);
        SumR_NF_RIS=sum(R_NF_RIS);
        SumR_FF_AP=sum(R_FF_AP);
        SumR_NF_AP=sum(R_NF_AP);
        SumR_MM=sum(R_MM);
        SumR_Rand=sum(R_Rand);
        
        MinR_FF_RIS=min(R_FF_RIS);
        MinR_NF_RIS=min(R_NF_RIS);
        MinR_FF_AP=min(R_FF_AP);
        MinR_NF_AP=min(R_NF_AP);
        MinR_MM=min(R_MM);
        MinR_Rand=min(R_Rand);
        
        temp_SumR_FF_RIS=temp_SumR_FF_RIS+SumR_FF_RIS./ITER;
        temp_SumR_NF_RIS=temp_SumR_NF_RIS+SumR_NF_RIS./ITER;
        temp_SumR_FF_AP=temp_SumR_FF_AP+SumR_FF_AP./ITER;
        temp_SumR_NF_AP=temp_SumR_NF_AP+SumR_NF_AP./ITER;
        temp_SumR_MM=temp_SumR_MM+SumR_MM./ITER;
        temp_SumR_Rand=temp_SumR_Rand+SumR_Rand./ITER;
        
        temp_MinR_FF_RIS=temp_MinR_FF_RIS+MinR_FF_RIS./ITER;
        temp_MinR_NF_RIS=temp_MinR_NF_RIS+MinR_NF_RIS./ITER;
        temp_MinR_FF_AP=temp_MinR_FF_AP+MinR_FF_AP./ITER;
        temp_MinR_NF_AP=temp_MinR_NF_AP+MinR_NF_AP./ITER;
        temp_MinR_MM=temp_MinR_MM+MinR_MM./ITER;
        temp_MinR_Rand=temp_MinR_Rand+MinR_Rand./ITER;
        
    end
    %% Lưu giá trị tốc độ sau khi kết thúc vòng lặp
    Record_SumR_FF_RIS(idx_K)=temp_SumR_FF_RIS;
    Record_SumR_NF_RIS(idx_K)=temp_SumR_NF_RIS;
    Record_SumR_FF_AP(idx_K)=temp_SumR_FF_AP;
    Record_SumR_NF_AP(idx_K)=temp_SumR_NF_AP;
    Record_SumR_MM(idx_K)=temp_SumR_MM;
    Record_SumR_Rand(idx_K)=temp_SumR_Rand;
    
    Record_MinR_FF_RIS(idx_K)=temp_MinR_FF_RIS;
    Record_MinR_NF_RIS(idx_K)=temp_MinR_NF_RIS;
    Record_MinR_FF_AP(idx_K)=temp_MinR_FF_AP;
    Record_MinR_NF_AP(idx_K)=temp_MinR_NF_AP;
    Record_MinR_MM(idx_K)=temp_MinR_MM;
    Record_MinR_Rand(idx_K)=temp_MinR_Rand;
    
end

xx=K_list;

%% Lưu file thông số
filename   =   strcat('xx',   '.mat');
save(['./',filename],    'xx','-v7.3');

filename   =   strcat('Record_SumR_FF_RIS',   '.mat');
save(['./',filename],    'Record_SumR_FF_RIS','-v7.3');

filename   =   strcat('Record_SumR_NF_RIS',   '.mat');
save(['./',filename],    'Record_SumR_NF_RIS','-v7.3');

filename   =   strcat('Record_SumR_FF_AP',   '.mat');
save(['./',filename],    'Record_SumR_FF_AP','-v7.3');

filename   =   strcat('Record_SumR_NF_AP',   '.mat');
save(['./',filename],    'Record_SumR_NF_AP','-v7.3');

filename   =   strcat('Record_SumR_MM',   '.mat');
save(['./',filename],    'Record_SumR_MM','-v7.3');

filename   =   strcat('Record_SumR_Rand',   '.mat');
save(['./',filename],    'Record_SumR_Rand','-v7.3');
%----------------------------------------------------%
filename   =   strcat('Record_MinR_FF_RIS',   '.mat');
save(['./',filename],    'Record_MinR_FF_RIS','-v7.3');

filename   =   strcat('Record_MinR_NF_RIS',   '.mat');
save(['./',filename],    'Record_MinR_NF_RIS','-v7.3');

filename   =   strcat('Record_MinR_FF_AP',   '.mat');
save(['./',filename],    'Record_MinR_FF_AP','-v7.3');

filename   =   strcat('Record_MinR_NF_AP',   '.mat');
save(['./',filename],    'Record_MinR_NF_AP','-v7.3');

filename   =   strcat('Record_MinR_MM',   '.mat');
save(['./',filename],    'Record_MinR_MM','-v7.3');

filename   =   strcat('Record_MinR_Rand',   '.mat');
save(['./',filename],    'Record_MinR_Rand','-v7.3');
%% Mì ăn liền,
%Nếu muốn xuất ra đồ thị liền thì bỏ comment rồi chạy section này tới hết code
%Từ đây tới hết code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Just Remove the comment of the following  lines, Run from this line to the end of this code
%load('xx.mat');
%load('Record_MinR_FF_RIS.mat');
%load('Record_MinR_NF_RIS.mat');
%load('Record_MinR_FF_AP.mat');
%load('Record_MinR_NF_AP.mat');
%load('Record_MinR_MM.mat');
%load('Record_MinR_Rand.mat');
% %-------------------------------%
%load('Record_SumR_FF_RIS.mat');
%load('Record_SumR_NF_RIS.mat');
%load('Record_SumR_FF_AP.mat');
%load('Record_SumR_NF_AP.mat');
%load('Record_SumR_MM.mat');
%load('Record_SumR_Rand.mat');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




label_fontsize=13;
legend_fontsize=11;
linewidth=1.5;
%% Đồ thị cho tốc độ tối thiểu
figure; hold on; box on; grid on;
p1 = plot(xx,Record_MinR_MM, 'Color', [246,83,20]./255,  'Linestyle', '-',  'Marker', '^'   , 'Linewidth', linewidth);
p2 = plot(xx,Record_MinR_FF_RIS, 'Color', [124,187,0]./255,  'Linestyle', '-',  'Marker', 's' , 'Linewidth', linewidth);
p3 = plot(xx,Record_MinR_NF_RIS, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', 'd'  , 'Linewidth', linewidth);
% p4 = plot(xx,Record_MinR_FF_AP, 'Color', [124,187,0]./255, 'Linestyle', '-',  'Marker', '<'    , 'Linewidth', linewidth);
% p5 = plot(xx,Record_MinR_NF_AP, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', '>'   , 'Linewidth', linewidth);
p6 = plot(xx,Record_MinR_Rand, 'Color', 255*[0.49,0.18,0.56]./255,  'Linestyle', '--',  'Marker', 'x'   , 'Linewidth', linewidth);

%set(gca,'XDir','reverse')
legend([p1,p3,p2,p6], {'Đề xuất đa tia', 'Thiết kế kênh tia cận trường','Thiết kế kênh tia viễn trường','Thiết kế pha ngẫu nhiên'},'fontsize',legend_fontsize)
xlabel('Số lượng người dùng','fontsize',label_fontsize)
ylabel('Tốc độ truyền tối thiểu[bits/s/Hz]','fontsize',label_fontsize)

toc


% %% Đồ thị tổng tốc độ
% figure; hold on; box on; grid on;
% p1 = plot(xx,Record_SumR_MM, 'Color', [246,83,20]./255,  'Linestyle', '-',  'Marker', '^'   , 'Linewidth', linewidth);
% p2 = plot(xx,Record_SumR_FF_RIS, 'Color', [124,187,0]./255,  'Linestyle', '-',  'Marker', 's' , 'Linewidth', linewidth);
% p3 = plot(xx,Record_SumR_NF_RIS, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', 'd'  , 'Linewidth', linewidth);
% p4 = plot(xx,Record_SumR_FF_AP, 'Color', [124,187,0]./255, 'Linestyle', '-',  'Marker', '<'    , 'Linewidth', linewidth);
% p5 = plot(xx,Record_SumR_NF_AP, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', '>'   , 'Linewidth', linewidth);
% p6 = plot(xx,Record_SumR_Rand, 'Color', 255*[0.49,0.18,0.56]./255,  'Linestyle', '--',  'Marker', 'x'   , 'Linewidth', linewidth);
% 
% set(gca,'XDir','reverse')
% legend([p1,p3,p2,p6], {'Proposed beam design', 'NFCB-based beam design','FFCB-based beam design','Random phase design'},'fontsize',legend_fontsize)
% xlabel('Number of UEs','fontsize',label_fontsize)
% ylabel('Achievable sum-rate [bits/s/Hz]','fontsize',label_fontsize)



