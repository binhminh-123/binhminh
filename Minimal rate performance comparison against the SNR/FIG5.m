clc;
clear all
close all
tic

%最小用户速率随SNR增长的变化

%% Para. Setup
N1 =64;
N2 = 8;
N = N1*N2; % the number of RIS elements
d = 0.5; % the antenna spacing 

%K是UE
K=8;
num_K=K;

% Iter Num. Fig in the paper is gene. with ITER = 600;
ITER = 100;

A = 4;%điều chỉnh Delta nhằm thay đổi kích thước codebook, Delta_1 càng lớn thì codebook càng ít mã
delta = 0.25;%?
D_oversample=1;%bước nhảy lấy mẫu

realsnr_max=0;
realsnr_list=[-10:1:realsnr_max];%mảng giá trị SNR


%Need Revise with different x
Bigfor_list=realsnr_list;
%% Gene Codebook
disp("Gene Near and Far Codebooks……")
% generate the far-field codebook
UN1 = exp(1i*2*pi*[0:(N1-1)*D_oversample]'*d*[0:N1-1]*(2/N1/D_oversample))/sqrt(N1);%ma trận số phức với phần thực là biên độ tín hiệu, phần ảo là pha của tín hiệu theo chiều ngang
UN2 = exp(1i*2*pi*[0:(N2-1)*D_oversample]'*d*[0:N2-1]*(2/N2/D_oversample))/sqrt(N2);%ma trận số phức với phần thực là biên độ tín hiệu, phần ảo là pha của tín hiệu theo chiều dọc
%do ko dùng góc trực tiếp nên dùng công thức [0:N1-1]*(2/N1/D_oversample) để tính góc chạy của tín hiệu và triệt tiêu lamda do = 1)
far_codebook = kron(UN1,UN2); %codebook là 1 mảng chứa các vector đại diện cho các hướng và cách thức truyền tín hiệu không dây

P3 = [2500*d,-2500*d,1200*d,200*d,0*d,-1000*d];%giới hạn vùng hoạt động theo kích thước của vùng Rayleigh
P4=P3;
P1=P3;
P2=P3;
Delta = 1*[100*d,100*d,100*d,100*d,100*d,100*d];%là vector chứa các giá trị khoảng cách
Delta1 = Delta*A;%dùng để tăng hoặc giảm kích thước codebook, tỷ lệ nghịch với A, khoảng mẫu
[near_codebook1,record] = generate_near_field_codebook(N1,N2,d,P3,P4,Delta1);%tạo codebook cho vùng cận trường
near_codebook1=near_codebook1./sqrt(N);%chuẩn hóa các vector trong codebook nhằm đảm bảo công suất truyền
disp("Finish Codebooks Gene")

%% Record the Rate(Final Used for Plot)
Record_SumR_FF_RIS=zeros(length(Bigfor_list),1);
Record_SumR_NF_RIS=zeros(length(Bigfor_list),1);
Record_SumR_FF_AP=zeros(length(Bigfor_list),1);
Record_SumR_NF_AP=zeros(length(Bigfor_list),1);
Record_SumR_MM=zeros(length(Bigfor_list),1);


Record_MinR_FF_RIS=zeros(length(Bigfor_list),1);
Record_MinR_NF_RIS=zeros(length(Bigfor_list),1);
Record_MinR_FF_AP=zeros(length(Bigfor_list),1);
Record_MinR_NF_AP=zeros(length(Bigfor_list),1);
Record_MinR_MM=zeros(length(Bigfor_list),1);


%% Main
t0 = clock;%bắt đầu tính thời gian chạy
for idx_Bigfor=1:length(Bigfor_list)
    LengthBigfor_list=length(Bigfor_list);
    %num_K=K_list(idx_K);
    
    %Need Revise with different x
    SNR_linear = 10.^(Bigfor_list(idx_Bigfor)/10.);%chuyển giá trị SNR từ giai dB sang giai tuyến tính để thực hiện tính toán
    %% Record buff (For the average calculate e.g. a=a+data./ITER)
    temp_SumR_FF_RIS=0;
    temp_SumR_NF_RIS=0;
    temp_SumR_FF_AP=0;
    temp_SumR_NF_AP=0;
    temp_SumR_MM=0;
   %tạo biến tạm để tí nữa tính giá trị trung bình
    
    temp_MinR_FF_RIS=0;
    temp_MinR_NF_RIS=0;
    temp_MinR_FF_AP=0;
    temp_MinR_NF_AP=0;
    temp_MinR_MM=0;
   
    %% Begin
    parfor idx_iter=1:ITER %chạy các vòng lặp song song nhau
        
        %% Gene Channel from BS to RIS
        %tạo mảng lưu các giá trị mã kênh
        FCCodewordsBuffer=zeros(N,num_K);
        NCCodewordsBuffer=zeros(N,num_K);
        PftCodewordsBuffer=zeros(N,num_K);
        %tạo mảng lưu các giá trị độ lợi kênh
        FCGainBuffer=zeros(num_K,1);
        NCGainBuffer=zeros(num_K,1);
        PftGainBuffer=zeros(num_K,1);
        %tạo mảng lưu giá trị mã lớp 2 có độ lợi lớn nhất
        max_index_2ndlayer=zeros(num_K,1);
        % generate the channel from the BS to the RIS
        %tạo kênh truyền từ BS tới RIS
        [G,px1,py1,pz1,alpha] = generate_G_near_field_channel(N1,N2,P1);
        GG=zeros(N,num_K);%tạo mảng lưu các giá trị kênh truyền cho từng người dùng (có 512 đường tới cho 1 người dùng)
        %% Record Time Clock
        %Need Revise with different x
        fprintf('For SNR  (NearField):i_num=%d of %d,i_iter=%d of %d | run %.4f s\n',idx_Bigfor,LengthBigfor_list,idx_iter,ITER,  etime(clock, t0));
        
        %% generate the channel from the RIS to the UE
        for k=1:num_K %tạo vòng lặp qua từng người dùng
            %tạo kênh truyền từ RIS-UE
            [hK,px2,py2,pz2,alpha] = generate_hr_near_field_channel(N1,N2,1,P2);%gọi hàm tạo kênh truyền cận trường từ RIS - UE
            %hK=hK./sqrt(N);%chuẩn hóa mảng vector hK
            Hc = diag(hK)*G;%tính các vector của mảng kênh truyền Hc bằng cách nhân ma trận chéo của hK với ma trận cột G, việc chuyển hK thành ma trận đường chéo giúp tính toán nhanh hơn
            GG(:,k)=Hc;%lưu các giá trị vector của mảng kênh Hc vào các cột ứng với từng người dùng từ UE1 - UE8
           
            %Far-Field BT%tìm ra kênh truyền có độ lợi lớn nhất cho từng người dùng
            [maxGainFC,idxFC]=max(abs(far_codebook*Hc).^2);%tìm giá trị độ lợi cao nhất trong mảng vào lưu vị trí của giá trị đó là idxFC
            FCCodewordsBuffer(:,k)=far_codebook(idxFC,:).';%lưu các từ mã(vector) tại hàng idxFC vào cột  k của mảng FCCodewordsBuffer
            FCGainBuffer(k)=maxGainFC;%lưu giá trị độ lợi lớn nhất vào mảng FCGainBuffer theo cột k
            
            %Near-Field BT
            array_gain = 0;%khởi tạo độ lợi mảng
            max_index=-1;
            for i =1:size(near_codebook1,1)%cho i chạy từ 1 đến hết số hàng của codebook cận trường
                if abs(near_codebook1(i,:)*Hc)^2>array_gain%so sánh độ lợi kênh truyền Hc với độ lợi mảng
                    max_index=i;%nếu lớn hơn thì cập nhật giá trị max_index = 1 
                    array_gain=abs(near_codebook1(i,:)*Hc)^2;%dùng giá trị độ lợi đó làm độ lới mới cho mảng
                end
            end
            NCCodewordsBuffer(:,k)=near_codebook1(max_index,:).';%lưu hàng có giá trị kênh độ lợi lớn nhất vào cột k của NCCodewordsBuffer
            % generate the second-level codes
            P21=[record(max_index,1)+Delta1(1)/2,record(max_index,1)-Delta1(1)/2,record(max_index,2)+Delta1(2)/2,record(max_index,2)-Delta1(2)/2,record(max_index,3)+Delta1(3)/2,record(max_index,3)-Delta1(3)/2];
            P22=[record(max_index,4)+Delta1(4)/2,record(max_index,4)-Delta1(4)/2,record(max_index,5)+Delta1(5)/2,record(max_index,5)-Delta1(5)/2,record(max_index,6)+Delta1(6)/2,record(max_index,6)-Delta1(6)/2];
            %thu hẹp vùng phân bố của BS và UE
            
            near_codebook2 = generate_near_field_codebook(N1,N2,d,P21,P22,Delta1*delta);%gọi hàm tạo codebook cận trường mới theo P21, P22 và Delta mới
            near_codebook2=near_codebook2./sqrt(N);%chuẩn hóa các vector trong codebook
            
            for i =1:size(near_codebook2,1)%cho i chạy hết codebook2
                if abs(near_codebook2(i,:)*Hc)^2>array_gain%so sánh độ lợi kênh cận trường mới với độ lợi mảng đã cập nhật lại ở trên
                    array_gain=abs(near_codebook2(i,:)*Hc)^2;%nếu lớn hơn thì gán làm độ lợi mảng mới
                    max_index_2ndlayer(k)=i;%lưu vị trí kênh lớp 2 đó
                end
            end
            %kiểm tra mã lớp 2 xem có hiệu quả hơn mã lớp 1
            if max_index_2ndlayer(k)>0 %nếu vị trí đó tồn tại mã lớp 2 hiệu quả thì cập nhật lại các giá trị độ lợi kênh và độ lợi mảng theo mã lớp 2
                NCCodewordsBuffer(:,k)=near_codebook2(max_index_2ndlayer(k),:).';
                maxGainNC=array_gain;
                NCGainBuffer(k)=abs(near_codebook2(max_index_2ndlayer(k),:)*Hc)^2;
            else
                %nếu không thì vẫn dùng giá trị theo codebook ban đầu
                NCGainBuffer(k)=abs(NCCodewordsBuffer(:,k).'*Hc)^2;
                
            end%
            
            %Tối ưu thành phần pha tín hiệu
            wc_opt = exp(1j*phase(Hc'));%tìm vector pha bằng cách phức hóa giá trị pha của kênh truyền Hc
            wc_opt=wc_opt./abs(wc_opt)/sqrt(N);%
            array_gainpft = abs(wc_opt*Hc)^2;%tính độ lợi mảng kênh truyền theo pha
            PftCodewordsBuffer(:,k)=wc_opt.';%ưu giá trị vectorpha vào cột k của mảng PftCodewordsBuffer
            PftGainBuffer(k)=array_gainpft;%lưu lại các giá trị độ lợi sau khi tối ưu pha tín hiệu
        end
        %% Deal with channel gain
        Product_mxg_DFT=prod(sqrt(FCGainBuffer));%tính tích độ lợi viễn trường của cả 8 người dùng(số thực)
        MultiBeamFC_Orig=FCCodewordsBuffer*((Product_mxg_DFT./sqrt(FCGainBuffer)));%tạo vector định hướng chùm tia viễn trường bằng cách nhân codeword viễn trường cho độ lợi viễn trường chung đã chuẩn hóa
        Product_mxg_NC=prod(sqrt(NCGainBuffer));%tính tích độ lợi cận trường của cả 8 người dùng
        MultiBeamNC_Orig=NCCodewordsBuffer*((Product_mxg_NC./sqrt(NCGainBuffer)));%tạo vector định hướng chùm tia cận trường bằng cách nhân codeword cận trường cho độ lợi cận trường chung đã chuẩn hóa
        Product_mxg_Pft=prod(sqrt(PftGainBuffer));%tính tích độ lợi cho cả 8 UE sau khi tối ưu hệ số pha
        MultiBeamPft_Orig=PftCodewordsBuffer*((Product_mxg_Pft./sqrt(PftGainBuffer)));%tạo vector định hướng cho chùm tia sau khi tối ưu hệ số pha
        %% Gene the Superpose FF-BF
        
        record_zeroFC=find(MultiBeamFC_Orig==0);%tìm các giá trị kênh bị trống
        MultiBeamFC_Orig(record_zeroFC)=exp(1j*2*pi*rand)/sqrt(N);%lấp các giá trị kênh trống đó với hàm ngẫu nhiên để đảm bảo không có tia nào bị rỗng(mảng vector)
        MultiBeamFCRIS=MultiBeamFC_Orig./abs(MultiBeamFC_Orig)/sqrt(N);%tạo đường đa tia tới mảng RIS và chuẩn hóa độ dài thành 1(mảng vector)
        MultiBeamFCAP=MultiBeamFC_Orig./max(abs(MultiBeamFC_Orig))/sqrt(N);%tạo đường đa tia tới UE, chọn ra đường tối ưu nhất và chuẩn hóa thành 1(mảng vector)
        %% Gene the Superpose NF-BF

        record_zeroNC=find(MultiBeamNC_Orig==0);
        MultiBeamNC_Orig(record_zeroNC)=exp(1j*2*pi*rand)/sqrt(N);(trong mảng MultiBeamNC_Orig, tại vị trí record_zeroNC, là vector exp(1j*2*pi*rand)/sqrt(N) ngẫu nhiên)
        MultiBeamNCRIS=MultiBeamNC_Orig./abs(MultiBeamNC_Orig)/sqrt(N);(mảng vector)
        MultiBeamNCDig=MultiBeamNC_Orig/norm(MultiBeamNC_Orig);%Use for MM's v_abs setting
        MultiBeamNCAP=MultiBeamNC_Orig./max(abs(MultiBeamNC_Orig))/sqrt(N);%Which means that the Amplitude & Phase can be adjusted. (But the Amplitude is in [0,1])
        %
        %% Gene the Superpose Pft-BF
        %MultiBeamPft_Orig=sum(PftCodewordsBuffer,2);
        record_zeroPft=find(MultiBeamPft_Orig==0);
        MultiBeamPft_Orig(record_zeroPft)=exp(1j*2*pi*rand)/sqrt(N);
        MultiBeamPftRIS=MultiBeamPft_Orig./abs(MultiBeamPft_Orig)/sqrt(N);
        MultiBeamPftDig=MultiBeamPft_Orig/norm(MultiBeamPft_Orig);%Use for MM's v_abs setting
        
        %
        %% Calculate the MultiBeam Gain
        G_MultiBeam_FFSuperpose=abs(MultiBeamFCRIS.'*GG).^2;%size:1 * K (RIS BF)
        G_MultiBeam_NFSuperpose=abs(MultiBeamNCRIS.'*GG).^2;%(RIS BF)
        G_MultiBeam_PftSuperpose=abs(MultiBeamPftDig.'*GG).^2;%(Dig BF)
        G_MultiBeam_FFSuperpose_AP=abs(MultiBeamFCAP.'*GG).^2;%size:1 * K (AP BF)
        G_MultiBeam_NFSuperpose_AP=abs(MultiBeamNCAP.'*GG).^2;%(AP BF)
        
        %% MM Algorithm
        loss3 = [];%khởi tạo mảng lưu giá trị suy hao
        eta=1;%hệ số kiểm soát độ lớn v_abs3
        v_abs3=max(sqrt(G_MultiBeam_PftSuperpose))*ones(1,num_K)'*eta;%vector chứa giá trị độ lớn của các phần tử trong vector v3(ngẫu nhiên)
        w3 = exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N);%là vector phức khởi tạo ngẫu nhiên và đc chuẩn hóa rồi cập nhật qua từng vòng lặp nhằm giảm thiểu khác biệt giữa v3 và 
        
        v_phase=exp(1j*1*pi*(2*rand(num_K,1)-1));%vector chứa pha dùng cho vector v3(ngẫu nhiên)
        v3 = v_abs3.*v_phase;%vector định hướng đc khởi tạo ngẫu nhiên để bắt đầu vòng lặp
        
        A=conj(NCCodewordsBuffer).';%chuyển vị liên hợp của từ mã cận trường
        A=diag(sqrt(NCGainBuffer))*conj(NCCodewordsBuffer).';%根据近场码本确定A     k*n
        %A=GG.';
        lambda = max(eig(A'*A));%tính giá trị riêng lớn nhất của tích ma trận A'.A
        max_iter=200;%giới hạn vòng lặp
        % Iteration
        for i=1:max_iter
            % update v_phase
            f = A*w3;%hàm thay thế
            v_phase = exp(1j * angle(f)); %cập nhật lại giá trị v_phase theo f
            v3 = v_abs3 .* v_phase;%cập nhật lại giá trị v3 theo v_phase mới
            % update w3
            w3_0iter=w3;
            for j_MM=1:10
                temp = A' * v3 - A' * A * w3_0iter + lambda * w3_0iter ;
                w3_0iter = exp( 1j*angle( temp ))/sqrt(N);
                %vòng lặp tối ưu giá trị w3
            end
            w3=w3_0iter;%gán giá trị w3 mới
            loss3 = [loss3, norm(v3 - A * w3, 2)];%cập nhật giá trị suy hao theo phép toán euclid bình phương cho thấy khác biệt giữa v3 và hàm f
        end
        %kết quả vòng lặp trả ra là vector hình thành tia w3 với giá trị biên độ v_abs3 và v_phase được tối ưu
        %% Calculate MM MultiBeam Gain
        G_MultiBeam_MM=abs(w3.'*GG).^2;
        G_MM=abs(w3.'*GG).^2;
        %% Calculate Rate for every UE
        R_FF_RIS=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose.');%Size K*1
        R_NF_RIS=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose.');
        R_FF_AP=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose_AP.');
        R_NF_AP=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose_AP.');
        R_MM=log2(1+SNR_linear.*G_MultiBeam_MM.');
        
        
        %% Calculate Sum-Rate & MinRate
        SumR_FF_RIS=sum(R_FF_RIS);
        SumR_NF_RIS=sum(R_NF_RIS);
        SumR_FF_AP=sum(R_FF_AP);
        SumR_NF_AP=sum(R_NF_AP);
        SumR_MM=sum(R_MM);
        
        
        MinR_FF_RIS=min(R_FF_RIS);
        MinR_NF_RIS=min(R_NF_RIS);
        MinR_FF_AP=min(R_FF_AP);
        MinR_NF_AP=min(R_NF_AP);
        MinR_MM=min(R_MM);
       

        
        temp_SumR_FF_RIS=temp_SumR_FF_RIS+SumR_FF_RIS./ITER;
        temp_SumR_NF_RIS=temp_SumR_NF_RIS+SumR_NF_RIS./ITER;
        temp_SumR_FF_AP=temp_SumR_FF_AP+SumR_FF_AP./ITER;
        temp_SumR_NF_AP=temp_SumR_NF_AP+SumR_NF_AP./ITER;
        temp_SumR_MM=temp_SumR_MM+SumR_MM./ITER;
        
        
        temp_MinR_FF_RIS=temp_MinR_FF_RIS+MinR_FF_RIS./ITER;
        temp_MinR_NF_RIS=temp_MinR_NF_RIS+MinR_NF_RIS./ITER;
        temp_MinR_FF_AP=temp_MinR_FF_AP+MinR_FF_AP./ITER;
        temp_MinR_NF_AP=temp_MinR_NF_AP+MinR_NF_AP./ITER;
        temp_MinR_MM=temp_MinR_MM+MinR_MM./ITER;
       
        
    end
    %% Record the Rate after Iter
    Record_SumR_FF_RIS(idx_Bigfor)=temp_SumR_FF_RIS;
    Record_SumR_NF_RIS(idx_Bigfor)=temp_SumR_NF_RIS;
    Record_SumR_FF_AP(idx_Bigfor)=temp_SumR_FF_AP;
    Record_SumR_NF_AP(idx_Bigfor)=temp_SumR_NF_AP;
    Record_SumR_MM(idx_Bigfor)=temp_SumR_MM;
    
    Record_MinR_FF_RIS(idx_Bigfor)=temp_MinR_FF_RIS;
    Record_MinR_NF_RIS(idx_Bigfor)=temp_MinR_NF_RIS;
    Record_MinR_FF_AP(idx_Bigfor)=temp_MinR_FF_AP;
    Record_MinR_NF_AP(idx_Bigfor)=temp_MinR_NF_AP;
    Record_MinR_MM(idx_Bigfor)=temp_MinR_MM;
    
            
       
end

xx=Bigfor_list;

%% Save
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



%% If you want to directly get the Fig,,use our local data.
%If you want to directly get the Fig., just Remove the comment and only run
%the following part from the next line to the end of this code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Just Remove the comment of the following  lines, Run from this line to the end of this code
% load('xx.mat');
% load('Record_MinR_FF_RIS.mat');
% load('Record_MinR_NF_RIS.mat');
% load('Record_MinR_FF_AP.mat');
% load('Record_MinR_NF_AP.mat');
% load('Record_MinR_MM.mat');
% 
% %-------------------------------%
% load('Record_SumR_FF_RIS.mat');
% load('Record_SumR_NF_RIS.mat');
% load('Record_SumR_FF_AP.mat');
% load('Record_SumR_NF_AP.mat');
% load('Record_SumR_MM.mat');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% END-If you want to directly get the Fig,use our local data.


label_fontsize=13;
legend_fontsize=11;
linewidth=1.5;

%% Figure for MinRate
figure; hold on; box on; grid on;
p1 = plot(xx,Record_MinR_MM, 'Color', [246,83,20]./255,  'Linestyle', '-',  'Marker', '^'   , 'Linewidth', linewidth);
p2 = plot(xx,Record_MinR_FF_RIS, 'Color', [124,187,0]./255,  'Linestyle', '-',  'Marker', 's' , 'Linewidth', linewidth);
p3 = plot(xx,Record_MinR_NF_RIS, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', 'd'  , 'Linewidth', linewidth);
% p4 = plot(xx,Record_MinR_FF_AP, 'Color', [124,187,0]./255, 'Linestyle', '-',  'Marker', '<'    , 'Linewidth', linewidth);
% p5 = plot(xx,Record_MinR_NF_AP, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', '>'   , 'Linewidth', linewidth);
% p6 = plot(xx,Record_MinR_Rand, 'Color', 255.*[0.49,0.18,0.56]./255,  'Linestyle', '--',  'Marker', 'x'  , 'Linewidth', linewidth);

%set(gca,'XDir','reverse')
legend([p1,p3,p2], {'Proposed beam design', 'NFCB-based beam design','FFCB-based beam design'},'fontsize',legend_fontsize)
xlabel('SNR [dB]','fontsize',label_fontsize)
ylabel('Minimal rate [bits/s/Hz]','fontsize',label_fontsize)

toc



