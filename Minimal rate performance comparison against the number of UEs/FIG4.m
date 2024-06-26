clc;
clear all
close all

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mô phỏng so sánh tỷ lệ truyền tối thiểu với số lượng người dùng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Thiết lập thông số mô phỏng
N1 =64; %số hàng của mảng phản xạ
N2 = 8; %số cột của mảng phản xạ
N = N1*N2; % số phần tử mảng RIS
d = 0.5; % khoảng cách giữa các phần tử: lambda/2


% K=4;
% U=K;
K_max=20;% số người dùng tối đa
K_list=[2:2:K_max];

% ITER là số lần lặp để tính trung bình, trong bài báo này thì ITER = 600;
ITER = 30; 

A = 4;
delta = 0.25;  %khoảng cách giữa các nút lưới codebook
D_oversample=1; %biến điều chỉnh khoảng phân bố các điểm lưới trong không gian của codebook

realsnr=5; %tỷ lệ tín hiệu đến nhiễu (SNR - Signal-to-Noise Ratio)
SNR_linear = 10.^(realsnr/10.); %SNR được chuyển đổi sang thang đo tuyến tính
%% Thiết lập codebook
disp("Gene Near and Far Codebooks……")
% Thiết lập codebook viễn trường
% Tính theo công thức 5a và 5b
UN1 = exp(1i*2*pi*[0:(N1-1)*D_oversample]'*d*[0:N1-1]*(2/N1/D_oversample))/sqrt(N1);
UN2 = exp(1i*2*pi*[0:(N2-1)*D_oversample]'*d*[0:N2-1]*(2/N2/D_oversample))/sqrt(N2);
far_codebook = kron(UN1,UN2); %Tích Kronecker giữa hai ma trận UN1⊗UN2

% Tọa độ phân bố của BS và user
P3 = [2500*d,-2500*d,1200*d,200*d,0*d,-1000*d];
P4=P3;
P1=P3;
P2=P3;
Delta = 1*[100*d,100*d,100*d,100*d,100*d,100*d]; %vector chứa các khoảng cách giữa người dùng và các phần tử trong mảng RIS
Delta1 = Delta*A; % nhân các giá trị vector Delta với giá trị A tạo ra vector Delta1 có thể điều chỉnh khoảng cách giữa các phần tử trong codebook cận trường.
[near_codebook1,record] = generate_near_field_codebook(N1,N2,d,P3,P4,Delta1);%near_codebook1: Mã codebook near field cho mảng RIS; record: ma trận tọa độ của các phần tử trong mảng RIS.
near_codebook1=near_codebook1./sqrt(N); %chuẩn hóa
disp("Finish Codebooks Gene")

%% Lưu số liệu tốc độ truyền ( để vẽ đồ thị )
%Hàm zeros(n, m) được sử dụng để tạo ma trận gồm n hàng và m cột với tất cả các phần tử được khởi tạo là 0
%length(K_list) cho biết số lượng phần tử của vector K_list, do đó số hàng sẽ bằng length(K_list), số cột là 1
Record_SumR_FF_RIS=zeros(length(K_list),1); %khởi tạo ma trận lưu tổng tốc độ truyền cho far field RIS
Record_SumR_NF_RIS=zeros(length(K_list),1); %khởi tạo ma trận lưu tổng tốc độ truyền cho near field RIS
Record_SumR_FF_AP=zeros(length(K_list),1); %khởi tạo ma trận lưu tổng tốc độ truyền cho mô hình AP far field
Record_SumR_NF_AP=zeros(length(K_list),1); %khởi tạo ma trận lưu tổng tốc độ truyền cho mô hình AP near field
Record_SumR_MM=zeros(length(K_list),1); %khởi tạo ma trận lưu tổng tốc độ truyền MM
Record_SumR_Rand=zeros(length(K_list),1); %khởi tạo ma trận lưu tổng tốc độ truyền cho thiết kế tia ngẫu nhiên

Record_MinR_FF_RIS=zeros(length(K_list),1); %khởi tạo ma trận lưu tốc độ truyền tối thiểu cho far field RIS
Record_MinR_NF_RIS=zeros(length(K_list),1); %khởi tạo ma trận lưu tốc độ truyền tối thiểu cho near field RIS
Record_MinR_FF_AP=zeros(length(K_list),1); %khởi tạo ma trận lưu tốc độ truyền tối thiểu cho mô hình AP far field
Record_MinR_NF_AP=zeros(length(K_list),1); %khởi tạo ma trận lưu tốc độ truyền tối thiểu cho mô hình AP near field
Record_MinR_MM=zeros(length(K_list),1); %khởi tạo ma trận lưu tốc độ truyền tối thiểu cho mô hình MM
Record_MinR_Rand=zeros(length(K_list),1); %khởi tạo ma trận lưu tốc độ truyền tối thiểu cho thiết kế tia ngẫu nhiên

%% Hàm chính
t0 = clock; %ghi lại thời điểm bắt đầu vòng lặp
for idx_K=1:length(K_list)
    LengthK_list=length(K_list);
    num_K=K_list(idx_K);
    K=num_K;
    %% Lưu biến tạm để tính trung bình (a=a+data./ITER), các biến tạm được khởi tạo với giá trị = 0
    temp_SumR_FF_RIS=0; %biến tạm để tính tổng tốc độ truyền của mô hình RIS với thiết kế tia viễn trường
    temp_SumR_NF_RIS=0; %biến tạm đểtính tổng tốc độ truyền của mô hình RIS với thiết kế tia cận trường
    temp_SumR_FF_AP=0; %biến tạm để tính tổng tốc độ truyền của mô hình AP với thiết kế tia viễn trường
    temp_SumR_NF_AP=0; %biến tạm để tính tổng tốc độ truyền của mô hình AP với thiết kế tia cận trường
    temp_SumR_MM=0; %biến tạm để tính tổng tốc độ truyền của mô hình MM
    temp_SumR_Rand=0; %biến tạm để tính tổng tốc độ truyền của mô hình với thiết kế tia ngẫu nhiên
    
    temp_MinR_FF_RIS=0; %biến tạm để tính tốc độ truyền tối thiểu của mô hình RIS với thiết kế tia viễn trường
    temp_MinR_NF_RIS=0; %biến tạm đểtính tốc độ truyền tối thiểu của mô hình RIS với thiết kế tia cận trường
    temp_MinR_FF_AP=0; %biến tạm để tính tốc độ truyền tối thiểu của mô hình AP với thiết kế tia viễn trường
    temp_MinR_NF_AP=0; %biến tạm để tính tốc độ truyền tối thiểu của mô hình AP với thiết kế tia cận trường
    temp_MinR_MM=0; %biến tạm để tính tốc độ truyền tối thiểu của mô hình MM
    temp_MinR_Rand=0; %biến tạm để tính tốc độ truyền tối thiểu của mô hình với thiết kế tia ngẫu nhiên
    %% Vòng lặp chính, mỗi lần lặp sẽ thực hiện một lần mô phỏng của kênh truyền từ BS tới RIS cho mỗi người dùng
    for idx_iter=1:ITER
        
        %% khởi tạo các kênh từ BS đến RIS
        FCCodewordsBuffer=zeros(N,num_K); %lưu mã truyền viễn trường
        NCCodewordsBuffer=zeros(N,num_K); %lưu mã truyền cận trường
        PftCodewordsBuffer=zeros(N,num_K); %lưu mã truyền tia precoding
        % lưu trữ các hệ số tín hiệu-tới-nhiễu (SNR) cho từng người dùng.
        FCGainBuffer=zeros(num_K,1); %vector lưu độ lợi viễn trường
        NCGainBuffer=zeros(num_K,1); %vector lưu độ lợi cận trường
        PftGainBuffer=zeros(num_K,1); %vector lưu độ lợi kênh precoding
        max_index_2ndlayer=zeros(num_K,1);
        % Tạo kênh truyền BS tới RIS
        [G,px1,py1,pz1,alpha] = generate_G_near_field_channel(N1,N2,P1);
        GG=zeros(N,num_K);
        %% Lưu thông số thời gian
        %số lần lặp hiện tại (idx_K), tổng số lần lặp (LengthK_list), số lần lặp hiện tại trong ITER (idx_iter), tổng số lần lặp trong ITER (ITER), và thời gian đã trôi qua kể từ thời điểm bắt đầu lặp (etime(clock, t0)
        fprintf('For UserNum  (NearField):i_num=%d of %d,i_iter=%d of %d | run %.4f s\n',idx_K,LengthK_list,idx_iter,ITER,  etime(clock, t0));
        
        %% Tạo kênh truyền RIS tới người dùng
        for k=1:num_K
            %gọi hàm để tạo kênh truyền từ RIS đến người dùng k và trả về ma trận kênh hK cùng với các thông số khác
            [hK,px2,py2,pz2,alpha] = generate_hr_near_field_channel(N1,N2,1,P2);
            %hK=hK./sqrt(N);
            Hc = diag(hK)*G; %được tính bằng cách nhân đường chéo của ma trận hK với ma trận G, và kết quả này được gán cho cột thứ k của ma trận GG.
            GG(:,k)=Hc;           
           %% Tạo tia cận trường và viễn trường：(Get N(or F)CCodewordsBuffer、N(or F)CGainBuffer、PftGainBuffer)
            % Tạo tia viễn trường
            [maxGainFC,idxFC]=max(abs(far_codebook*Hc).^2); %tính giá trị lớn nhất của module bình phương kênh truyền từ RIS tới người dùng hiện tại
            FCCodewordsBuffer(:,k)=far_codebook(idxFC,:).'; %lưu giá trị trên vào cột thứ k của ma trận FCCodewordsBuffer: vector precoding được tạo ra cho kênh truyền BS đến RIS
            FCGainBuffer(k)=maxGainFC; %biến lưu giá trị maxgainFC cho người dùng thứ k
            
            % Tạo tia cận trường
            array_gain = 0; %Biến này sẽ lưu trữ giá trị lớn nhất của module bình phương kênh truyền từ RIS đến người dùng.
            max_index=-1; %Biến này sẽ lưu trữ chỉ số của vectơ codebook trong codebook cận trường với người dùng.
            for i =1:size(near_codebook1,1)
                if abs(near_codebook1(i,:)*Hc)^2>array_gain
                    max_index=i;
                    array_gain=abs(near_codebook1(i,:)*Hc)^2;
                end
            end
            NCCodewordsBuffer(:,k)=near_codebook1(max_index,:).'; %Gán vectơ codebook tìm được vào NCCodewordsBuffer (vector precoding được tạo ra cho kênh truyền từ RIS đến người dùng)
            % tạo ra tọa độ vị trí của người dùng hiện tại dựa trên vectơ codebook tối ưu trong near_codebook1. Tọa độ được tính dựa trên tọa độ trung tâm và độ lệch từ vectơ codebook được chọn và các giá trị trong vectơ Delta1.
            %Tạo ra các tọa độ theo chiều dọc của người dùng hiện tại
            P21=[record(max_index,1)+Delta1(1)/2,record(max_index,1)-Delta1(1)/2,record(max_index,2)+Delta1(2)/2,record(max_index,2)-Delta1(2)/2,record(max_index,3)+Delta1(3)/2,record(max_index,3)-Delta1(3)/2];
            %Tạo ra các tọa độ theo chiều ngang của người dùng hiện tại
            P22=[record(max_index,4)+Delta1(4)/2,record(max_index,4)-Delta1(4)/2,record(max_index,5)+Delta1(5)/2,record(max_index,5)-Delta1(5)/2,record(max_index,6)+Delta1(6)/2,record(max_index,6)-Delta1(6)/2];

            %tạo ra 1 codebook cận trường thứ 2
            near_codebook2 = generate_near_field_codebook(N1,N2,d,P21,P22,Delta1*delta);
            near_codebook2=near_codebook2./sqrt(N); %chia cho căn N để chuẩn hóa        
            for i =1:size(near_codebook2,1)
                if abs(near_codebook2(i,:)*Hc)^2>array_gain
                    array_gain=abs(near_codebook2(i,:)*Hc)^2;
                    max_index_2ndlayer(k)=i;
                end
            end
            if max_index_2ndlayer(k)>0 %kiểm tra nếu chỉ số của hàng có giá trị lớn nhất trong near_codebook2 lớn hơn 0. Nếu đúng, có nghĩa là cận trường thứ hai đã được áp dụng
                NCCodewordsBuffer(:,k)=near_codebook2(max_index_2ndlayer(k),:).';
                maxGainNC=array_gain;
                NCGainBuffer(k)=abs(near_codebook2(max_index_2ndlayer(k),:)*Hc)^2;
            else        
                NCGainBuffer(k)=abs(NCCodewordsBuffer(:,k).'*Hc)^2; %biến cho biết mức độ mạnh yếu của tín hiệu mà người dùng nhận được
            end
            
            %Pft precoding
            %%dòng dưới tạo ra ma trận omega theo công thức (21)
            wc_opt = exp(1j*phase(Hc')); %If can not work in this line maybe with different Matlab version, you can try to use "angle" function to replace "phase"
            wc_opt=wc_opt./abs(wc_opt)/sqrt(N); %chuẩn hóa ma trận trên
            array_gainpft = abs(wc_opt*Hc)^2; %cho biết mức độ mạnh yếu của tín hiệu truyền qua kênh từ RIS đến người dùng sau khi áp dụng PFT (Precoding)
            PftCodewordsBuffer(:,k)=wc_opt.'; %gán giá trị wc_opt đã được chuẩn hóa vào cột thứ k của ma trận PftCodewordsBuffer (dùng kỹ thuật PFT để để tăng cường hiệu suất truyền thông từ BS đến RIS hoặc từ RIS đến người dùng)
            PftGainBuffer(k)=array_gainpft; %lưu giá trị sau khi dùng PFT cho người dùng thứ k vào PftGainBuffer.
        end
        %% Xử lý độ lợi kênh truyền: tính các vector để dùng cho precoding
        Product_mxg_DFT=prod(sqrt(FCGainBuffer)); %hàm prod dùng để tính tích của căn bậc 2 các phần tử của ma trận
        MultiBeamFC_Orig=FCCodewordsBuffer*((Product_mxg_DFT./sqrt(FCGainBuffer))); %biến để lưu tín hiệu sau khi đã qua precoding
        Product_mxg_NC=prod(sqrt(NCGainBuffer));
        MultiBeamNC_Orig=NCCodewordsBuffer*((Product_mxg_NC./sqrt(NCGainBuffer)));
        Product_mxg_Pft=prod(sqrt(PftGainBuffer));
        MultiBeamPft_Orig=PftCodewordsBuffer*((Product_mxg_Pft./sqrt(PftGainBuffer)));
        %% Chồng chập tia viễn trường
        %MultiBeamFC_Orig=sum(FCCodewordsBuffer,2);
        record_zeroFC=find(MultiBeamFC_Orig==0);%tìm các vị trí trong biến MultiBeamFC_Orig mà giá trị là 0.
        MultiBeamFC_Orig(record_zeroFC)=exp(1j*2*pi*rand)/sqrt(N);%thay thế các giá trị 0 tìm được bằng một giá trị ngẫu nhiên được tạo ra từ phân phối đều trên đơn vị phức và được chuẩn hóa bởi căn bậc hai của N.
        MultiBeamFCRIS=MultiBeamFC_Orig./abs(MultiBeamFC_Orig)/sqrt(N);%chuẩn hóa các giá trị trong MultiBeamFC_Orig để tạo ra một tín hiệu với biên độ bằng 1 và được sử dụng cho RIS.
        MultiBeamFCAP=MultiBeamFC_Orig./max(abs(MultiBeamFC_Orig))/sqrt(N);% %chuẩn hóa các giá trị trong MultiBeamFC_Orig để tạo ra một tín hiệu với biên độ tối đa là 1 và được sử dụng cho AP (Antenna Precoding)
        %% Chồng chập tia cận trường (các biến được thực hiện tương tự như viễn trường)
        %MultiBeamNC_Orig=sum(NCCodewordsBuffer,2);
        record_zeroNC=find(MultiBeamNC_Orig==0);
        MultiBeamNC_Orig(record_zeroNC)=exp(1j*2*pi*rand)/sqrt(N);
        MultiBeamNCRIS=MultiBeamNC_Orig./abs(MultiBeamNC_Orig)/sqrt(N);
        MultiBeamNCDig=MultiBeamNC_Orig/norm(MultiBeamNC_Orig);%Use for MM's v_abs setting
        MultiBeamNCAP=MultiBeamNC_Orig./max(abs(MultiBeamNC_Orig))/sqrt(N);%Which means that the Amplitude & Phase can be adjusted. (But the Amplitude is in [0,1])
        %
        %% Chồng chập tia Pft ( kết hợp Tiền xử lý và Kết hợp tín hiệu vào tạo tia )
        %MultiBeamPft_Orig=sum(PftCodewordsBuffer,2);
        record_zeroPft=find(MultiBeamPft_Orig==0); %Tìm các vị trí trong MultiBeamPft_Orig mà giá trị = 0 
        MultiBeamPft_Orig(record_zeroPft)=exp(1j*2*pi*rand)/sqrt(N); %thay các vị trí tìm được bằng các số phức ngẫu nhiên được tạo ra từ phân phối đều trên khoảng từ 0 đến 2π.
        MultiBeamPftRIS=MultiBeamPft_Orig./abs(MultiBeamPft_Orig)/sqrt(N); %Tính toán vectơ các tia sau khi chồng chập và chi cho abs của nó để chuẩn hóa
        MultiBeamPftDig=MultiBeamPft_Orig/norm(MultiBeamPft_Orig);%Use for MM's v_abs setting, vẫn là vecto ở trên nhưng chuẩn hóa theo độ dài của nó (norm)
        MultiBeamRand=exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N); %Tạo vectơ tia phát ngẫu nhiên với phân phối đều và chuẩn hóa nó.
        %
        %% Tính độ lợi đa tia
        G_MultiBeam_FFSuperpose=abs(MultiBeamFCRIS.'*GG).^2;%size:1 * K (RIS BF), độ lợi đa tia cho far field khi sử dụng kỹ thuật beamforming dựa trên tín hiệu từ RIS
        G_MultiBeam_NFSuperpose=abs(MultiBeamNCRIS.'*GG).^2;%(RIS BF), như trên nhưng cho near field
        G_MultiBeam_PftSuperpose=abs(MultiBeamPftDig.'*GG).^2;%(Dig BF), sử dụng thêm phương pháp pft vào
        G_MultiBeam_FFSuperpose_AP=abs(MultiBeamFCAP.'*GG).^2;%size:1 * K (AP BF), như trên nhưng dựa vào tín hiệu từ RIS và AP
        G_MultiBeam_NFSuperpose_AP=abs(MultiBeamNCAP.'*GG).^2;%(AP BF), như trên nhưng cho near field
        G_MultiBeam_Rand=abs(MultiBeamRand.'*GG).^2;%(AP BF), dựa theo tín hiệu được phát sóng ngẫu nhiên

        %% Thuật toán MM
        loss3 = [];% mảng để lưu trữ giá trị của hàm suy hao trong quá trình tối ưu hoá.
        eta=1;
        v_abs3=max(sqrt(G_MultiBeam_PftSuperpose))*ones(1,num_K)'*eta; %vector amplitude của tín hiệu truyền đi v_abs3 dựa trên giá trị lớn nhất của độ lợi đa tia từ PFT và nhân với hệ số eta.
        w3 = exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N);%vector trọng số của RIS được khởi tạo ngẫu nhiên.
        
        v_phase=exp(1j*1*pi*(2*rand(num_K,1)-1));%vector pha của RIS được khởi tạo ngẫu nhiên
        v3 = v_abs3.*v_phase;
        
        A=conj(NCCodewordsBuffer).';% ma trận chuyển vị liên hợp của ma trận codebook cận trường
        A=diag(sqrt(NCGainBuffer))*conj(NCCodewordsBuffer).';% nhân ma trận trên với ma trận đường chéo có các phần tử là căn bậc 2 của NCGainBuffer
                
        %A=GG.';
        lambda = max(eig(A'*A));%tính giá trị lớn nhất của giá trị riêng của ma trận A'*A
        max_iter=200; %số lần lặp tối đa
        % Vòng lặp
        for i=1:max_iter
            % update v_phase
            f = A*w3;
            v_phase = exp(1j * angle(f)); % cập nhật biến phụ bằng cách lấy góc của phần thực của ma trận f.
            v3 = v_abs3 .* v_phase;
            % cập nhật thông số w3
            w3_0iter=w3;
            for j_MM=1:10
                temp = A' * v3 - A' * A * w3_0iter + lambda * w3_0iter ;% tính biến tạm
                w3_0iter = exp( 1j*angle( temp ))/sqrt(N);  
            end
            w3=w3_0iter;
            loss3 = [loss3, norm(v3 - A * w3, 2)];
        end
        %% Tính độ lợi đa tia bằng MM
        G_MultiBeam_MM=abs(w3.'*GG).^2; %công suất tín hiệu thu được tại các người dùng khi sử dụng phương pháp MM
        G_MM=abs(w3.'*GG).^2;%giống cái trên
        %% Tính tốc độ truyền cho từng người dùng
        R_FF_RIS=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose.');%Size K*1 % khi sử dụng beamforming RIS viễn trường
        R_NF_RIS=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose.');% khi sử dụng beamforming RIS cận trường
        R_FF_AP=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose_AP.');%khi sử dụng beamforming AP viễn trường
        R_NF_AP=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose_AP.');%khi sử dụng beamforming AP cận trường
        R_MM=log2(1+SNR_linear.*G_MultiBeam_MM.');% khi sử dụng thuật toán MM
        R_Rand=log2(1+SNR_linear.*G_MultiBeam_Rand.');% khi dùng độ lợi ngẫu nhiên
        %% Tổng tốc độ truyền 
        SumR_FF_RIS=sum(R_FF_RIS);
        SumR_NF_RIS=sum(R_NF_RIS);
        SumR_FF_AP=sum(R_FF_AP);
        SumR_NF_AP=sum(R_NF_AP);
        SumR_MM=sum(R_MM);
        SumR_Rand=sum(R_Rand);
        %Tốc độ truyền tối thiểu
        MinR_FF_RIS=min(R_FF_RIS);
        MinR_NF_RIS=min(R_NF_RIS);
        MinR_FF_AP=min(R_FF_AP);
        MinR_NF_AP=min(R_NF_AP);
        MinR_MM=min(R_MM);
        MinR_Rand=min(R_Rand);
        % tổng trung bình của tốc độ truyền cho từng trường hợp (FF_RIS, NF_RIS, FF_AP, NF_AP, MM, và Random) qua các lần lặp
        temp_SumR_FF_RIS=temp_SumR_FF_RIS+SumR_FF_RIS./ITER;
        temp_SumR_NF_RIS=temp_SumR_NF_RIS+SumR_NF_RIS./ITER;
        temp_SumR_FF_AP=temp_SumR_FF_AP+SumR_FF_AP./ITER;
        temp_SumR_NF_AP=temp_SumR_NF_AP+SumR_NF_AP./ITER;
        temp_SumR_MM=temp_SumR_MM+SumR_MM./ITER;
        temp_SumR_Rand=temp_SumR_Rand+SumR_Rand./ITER;
        % tính giá trị trung bình của tốc độ truyền tối thiểu cho từng trường hợp qua các lần lặp
        temp_MinR_FF_RIS=temp_MinR_FF_RIS+MinR_FF_RIS./ITER;
        temp_MinR_NF_RIS=temp_MinR_NF_RIS+MinR_NF_RIS./ITER;
        temp_MinR_FF_AP=temp_MinR_FF_AP+MinR_FF_AP./ITER;
        temp_MinR_NF_AP=temp_MinR_NF_AP+MinR_NF_AP./ITER;
        temp_MinR_MM=temp_MinR_MM+MinR_MM./ITER;
        temp_MinR_Rand=temp_MinR_Rand+MinR_Rand./ITER;
        
    end
    %% Lưu các giá trị tốc độ sau khi kết thúc vòng lặp
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
%Các lệnh save này tạo ra các file .mat chứa các biến và cấu trúc dữ liệu tương ứng, với định dạng nén -v7.3 để hỗ trợ việc lưu trữ dữ liệu lớn. 
%hàm strcat để nối các chuỗi lại với nhau
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

%Nếu muốn xuất ra đồ thị liền thì bỏ comment rồi chạy section này tới hết code
%Từ đây tới hết code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



