clc;
clear all
close all
tic %bắt đầu tính thời gian chạy của code

%Mô phỏng so sánh tốc độ truyền tối thiểu với SNR

%% Cấu hình thông số mô phỏng
N1 =64; %số hàng của mảng phản xạ
N2 = 8; %số cột của mảng phản xạ
N = N1*N2; % số phần tử RIS
d = 0.5; % khoảng cách giữa các phần tử lambda/2


K=8; %số người dùng cuối
num_K=K;

% ITER là số lần lặp để tính kết quả trung bình, trong bài báo này ITER = 600;
ITER = 30;

A = 4;
delta = 0.25; %khoảng cách giữa các nút lưới codebook
D_oversample=1; %biến điều chỉnh khoảng phân bố các điểm lưới trong không gian của codebook

realsnr_max=0; % giá trị nhiễu tối đa
realsnr_list=[-10:1:realsnr_max]; % khoảng giá trị SNR



Bigfor_list=realsnr_list; %gán giá trị SNR cho biến phụ Bigfor_list

%% Tạo codebook
disp("Gene Near and Far Codebooks……")

% Tạo codebook viễn trường
UN1 = exp(1i*2*pi*[0:(N1-1)*D_oversample]'*d*[0:N1-1]*(2/N1/D_oversample))/sqrt(N1); % tạo codebook viễn trường cho chiều thứ 1
UN2 = exp(1i*2*pi*[0:(N2-1)*D_oversample]'*d*[0:N2-1]*(2/N2/D_oversample))/sqrt(N2); % tạo codebook viễn trường cho chiều thứ 2
far_codebook = kron(UN1,UN2); % phép nhân Kronecker

P3 = [2500*d,-2500*d,1200*d,200*d,0*d,-1000*d]; % tạo vector P3 chứa các giá trị thực
P4=P3; % sao chép P3 vào biến phụ P4
P1=P3; % sao chép P3 vào biến phụ P1
P2=P3; % sao chép P3 vào biến phụ P2
Delta = 1*[100*d,100*d,100*d,100*d,100*d,100*d];
Delta1 = Delta*A; % nhân các giá trị vector Delta với giá trị A tạo ra vector Delta1 có thể điều chỉnh khoảng cách giữa các phần tử trong codebook cận trường.
[near_codebook1,record] = generate_near_field_codebook(N1,N2,d,P3,P4,Delta1); % tạo codebook cận trường
near_codebook1=near_codebook1./sqrt(N); % chuẩn hóa bằng chia căn bậc 2 N
disp("Finish Codebooks Gene")

%% Lưu thông số tốc độ truyền để vẽ đồ thị
Record_SumR_FF_RIS=zeros(length(Bigfor_list),1); %lưu giá trị trung bình cộng tốc độ truyền của mô hình RIS với thiết kế tia viễn trường
Record_SumR_NF_RIS=zeros(length(Bigfor_list),1); %lưu giá trị trung bình cộng tốc độ truyền của mô hình RIS với thiết kế tia cận trường
Record_SumR_FF_AP=zeros(length(Bigfor_list),1); %lưu giá trị trung bình cộng tốc độ truyền của mô hình AP với thiết kế tia viễn trường
Record_SumR_NF_AP=zeros(length(Bigfor_list),1); %lưu giá trị trung bình cộng tốc độ truyền của mô hình AP với thiết kế tia cận trường
Record_SumR_MM=zeros(length(Bigfor_list),1); %lưu giá trị trung bình cộng tốc độ truyền của mô hình MM
Record_SumR_Rand=zeros(length(Bigfor_list),1); %lưu giá trị trung bình cộng tốc độ truyền của mô hình với thiết kế tia ngẫu nhiên

Record_MinR_FF_RIS=zeros(length(Bigfor_list),1); %lưu tốc độ truyền tối thiểu của mô hình RIS với thiết kế tia viễn trường
Record_MinR_NF_RIS=zeros(length(Bigfor_list),1); %lưu tốc độ truyền tối thiểu của mô hình RIS với thiết kế tia cận trường
Record_MinR_FF_AP=zeros(length(Bigfor_list),1); %lưu tốc độ truyền tối thiểu của mô hình AP với thiết kế tia viễn trường
Record_MinR_NF_AP=zeros(length(Bigfor_list),1); %lưu tốc độ truyền tối thiểu của mô hình AP với thiết kế tia cận trường
Record_MinR_MM=zeros(length(Bigfor_list),1); %lưu tốc độ truyền tối thiểu của mô hình MM
Record_MinR_Rand=zeros(length(Bigfor_list),1); %lưu tốc độ truyền tối thiểu của mô hình với thiết kế tia ngẫu nhiên


%% Hàm chính
t0 = clock; % ghi lại thời gian bắt đầu vòng lặp
for idx_Bigfor=1:length(Bigfor_list) % bắt đầu vòng lặp qua các giá trị Bigfor_list. idx_Bigfor là chỉ số của phần tử hiện tại trong vòng lặp
    LengthBigfor_list=length(Bigfor_list); %gán độ dài của Bigfor_list vào biến LengthBigfor_list
    %num_K=K_list(idx_K);
    
    
    SNR_linear = 10.^(Bigfor_list(idx_Bigfor)/10.); %Thay đổi giá trị SNR từ giai dB sang tuyến tính

    %% Lưu bộ đệm để tính trung bình (a=a+data./ITER)
    temp_SumR_FF_RIS=0;%biến tạm để tính giá trị trung bình cộng tốc độ truyền của mô hình RIS với thiết kế tia viễn trường
    temp_SumR_NF_RIS=0;%biến tạm đểtính giá trị trung bình cộng tốc độ truyền của mô hình RIS với thiết kế tia cận trường
    temp_SumR_FF_AP=0;%biến tạm để tính giá trị trung bình cộng tốc độ truyền của mô hình AP với thiết kế tia viễn trường
    temp_SumR_NF_AP=0;%biến tạm để tính giá trị trung bình cộng tốc độ truyền của mô hình AP với thiết kế tia cận trường
    temp_SumR_MM=0;%biến tạm để tính giá trị trung bình cộng tốc độ truyền của mô hình MM
    temp_SumR_Rand=0;%biến tạm để tính giá trị trung bình cộng tốc độ truyền của mô hình với thiết kế tia ngẫu nhiên

    
    temp_MinR_FF_RIS=0;%biến tạm để tính tốc độ truyền tối thiểu của mô hình RIS với thiết kế tia viễn trường
    temp_MinR_NF_RIS=0;%biến tạm đểtính tốc độ truyền tối thiểu của mô hình RIS với thiết kế tia cận trường
    temp_MinR_FF_AP=0;%biến tạm để tính tốc độ truyền tối thiểu của mô hình AP với thiết kế tia viễn trường
    temp_MinR_NF_AP=0;%biến tạm để tính tốc độ truyền tối thiểu của mô hình AP với thiết kế tia cận trường
    temp_MinR_MM=0;%biến tạm để tính tốc độ truyền tối thiểu của mô hình MM
    temp_MinR_Rand=0;%biến tạm để tính tốc độ truyền tối thiểu của mô hình với thiết kế tia ngẫu nhiên

    %% Vòng lặp chính
    parfor idx_iter=1:ITER
        
        %% Tạo ma trận lưu các giá trị mã và độ lợi BS tới RIS
        FCCodewordsBuffer=zeros(N,num_K); % khởi tạo ma trận với kích thước N x num_K để lưu các mã truyền tia viễn trường
        NCCodewordsBuffer=zeros(N,num_K); % khởi tạo ma trận với kích thước N x num_K để lưu các mã truyền tia cận trường
        PftCodewordsBuffer=zeros(N,num_K); % khởi tạo ma trận với kích thước N x num_K để lưu các mã truyền tia kỹ thuật Precoding(PFT)
        FCGainBuffer=zeros(num_K,1); % khởi tạo vector với kich thước num_K x 1 để lưu độ lợi kênh viễn trường
        NCGainBuffer=zeros(num_K,1); % khởi tạo vector với kich thước num_K x 1 để lưu độ lợi kênh cận trường
        PftGainBuffer=zeros(num_K,1); % khởi tạo vector với kich thước num_K x 1 để lưu độ lợi kênh kỹ thuật Precoding(PFT)
        max_index_2ndlayer=zeros(num_K,1); % khởi tạo vector có kích thước num_K x 1 để lưu chỉ số tốt nhất của mã truyền tín hiệu trong kỹ thuật mã 2 lớp

        % tạo kênh BS - RIS
        [G,px1,py1,pz1,alpha] = generate_G_near_field_channel(N1,N2,P1); % tạo kênh truyền từ trạm cơ sở (BS) tới bề mặt phản xạ thông qua hệ thống Reflector Intelligent Surface (RIS) trong mô hình gần trường
        GG=zeros(N,num_K); % tạo vector để lưu kết quả cuối cùng của kênh truyền




        %% Lưu thông số thời gian
        %In ra các thông tin trạng thái của quá trình tính toán
        fprintf('For SNR  (NearField):i_num=%d of %d,i_iter=%d of %d | run %.4f s\n',idx_Bigfor,LengthBigfor_list,idx_iter,ITER,  etime(clock, t0));
        
        %% Tạo kênh truyền từ RIS tới UE
        for k=1:num_K
            [hK,px2,py2,pz2,alpha] = generate_hr_near_field_channel(N1,N2,1,P2); % gọi hàm tạo kênh cận trường
            %hK=hK./sqrt(N); % dùng để chuẩn hóa vector hK
            Hc = diag(hK)*G; % tính toán kênh truyền từ RIS đến người dùng cuối bằng cách nhân ma trận chéo với vector kênh hK và sau đó nhân kết quả với ma trận G
            GG(:,k)=Hc; % kênh truyền từ RIS đến người dùng cuối được lưu vào cột thứ k của ma trận GG
            %% FF and NF BT：(Get N(or F)CCodewordsBuffer、N(or F)CGainBuffer、PftGainBuffer)
            %Tạo tia viễn trường
            [maxGainFC,idxFC]=max(abs(far_codebook*Hc).^2); % tìm giá trị lớn nhất của bình phương của độ tuyệt đối tích vô hướng, đại diện cho độ lợi lớn nhất của các mã tia viễn trường
            FCCodewordsBuffer(:,k)=far_codebook(idxFC,:).'; % lưu trữ mã tia viễn trường tốt nhất vào FCCodewordsBuffer cho người dùng thứ k
            FCGainBuffer(k)=maxGainFC; % lưu trữ độ lợi lớn nhất vào FCGainBuffer cho người dùng thứ k
            
            %Tạo tia cận trường
            array_gain = 0; % khởi tạo ban đầu với giá trị 0 để lưu trữ độ lợi lớn nhất tìm thấy
            max_index=-1; % khởi tạo ban đầu với giá trị -1 để lưu trữ chỉ số của mã tia cận trường có độ lợi lớn nhất
            for i =1:size(near_codebook1,1) % Vòng lặp for duyệt qua tất cả các mã trong near_codebook1
                if abs(near_codebook1(i,:)*Hc)^2>array_gain % tính bình phương của giá trị tuyệt đối của tích vô hướng giữa mã cận trường near_codebook1(i,:) và kênh truyền từ RIS đến người dùng thứ k (Hc)
                    max_index=i; % Sau khi hoàn thành vòng lặp, max_index sẽ chứa chỉ số của mã tia cận trường có độ lợi lớn nhất
                    array_gain=abs(near_codebook1(i,:)*Hc)^2; % Nếu giá trị độ lợi mới tìm thấy lớn hơn array_gain, thì cập nhật array_gain và max_index
                end
            end
            NCCodewordsBuffer(:,k)=near_codebook1(max_index,:).'; % Mã tia cận trường tốt nhất được lưu vào NCCodewordsBuffer cho người dùng thứ k

            % Tạo mã 2 lớp
            P21=[record(max_index,1)+Delta1(1)/2,record(max_index,1)-Delta1(1)/2,record(max_index,2)+Delta1(2)/2,record(max_index,2)-Delta1(2)/2,record(max_index,3)+Delta1(3)/2,record(max_index,3)-Delta1(3)/2];
            P22=[record(max_index,4)+Delta1(4)/2,record(max_index,4)-Delta1(4)/2,record(max_index,5)+Delta1(5)/2,record(max_index,5)-Delta1(5)/2,record(max_index,6)+Delta1(6)/2,record(max_index,6)-Delta1(6)/2];
            % các mảng P21 và P22 chứa tọa độ của các điểm mẫu cho 2 lớp trong thiết kế mã tia cận trường, các điểm mẫu này được xác định bằng cách thêm và trừ một nửa giá trị của Delta1 từ và đến các tọa độ tương ứng của mã tia cận trường đã chọn (record(max_index,:))
             
            near_codebook2 = generate_near_field_codebook(N1,N2,d,P21,P22,Delta1*delta); % gọi hàm tạo codebook
            near_codebook2=near_codebook2./sqrt(N); % bước chuẩn hóa
            
            for i =1:size(near_codebook2,1) % lặp qua các hàng của ma trận near_codebook2 và tính toán độ lợi của mỗi vectơ cột trong ma trận khi nhân với vectơ kênh Hc
                if abs(near_codebook2(i,:)*Hc)^2>array_gain % Nếu độ lợi này lớn hơn array_gain, thì cập nhật array_gain và chỉ số max_index_2ndlayer(k) của hàng tương ứng
                    array_gain=abs(near_codebook2(i,:)*Hc)^2;
                    max_index_2ndlayer(k)=i;
                end
            end
            if max_index_2ndlayer(k)>0 % kiểm tra xem mã 2 lớp đã chạy chưa
                NCCodewordsBuffer(:,k)=near_codebook2(max_index_2ndlayer(k),:).'; % Nếu đã chạy, nó sẽ sử dụng vectơ mã từ near_codebook2 tại chỉ số max_index_2ndlayer(k) và tính toán độ lợi của vectơ này
                maxGainNC=array_gain;
                NCGainBuffer(k)=abs(near_codebook2(max_index_2ndlayer(k),:)*Hc)^2;
            else
                
                NCGainBuffer(k)=abs(NCCodewordsBuffer(:,k).'*Hc)^2; % Nếu chưa chạy, nó sẽ tính toán độ lợi của vectơ mã từ NCCodewordsBuffer tại cùng một chỉ số k
                
            end%
            
            %Kỹ thuật Precoding(PFT)
            wc_opt = exp(1j*phase(Hc')); % tính toán vectơ trọng số pha tối ưu (wc_opt) dựa trên pha của kênh (Hermitian conjugate của Hc)
            wc_opt=wc_opt./abs(wc_opt)/sqrt(N); % chuẩn hóa
            array_gainpft = abs(wc_opt*Hc)^2; %  tính toán độ lợi của tín hiệu sau khi áp dụng phương pháp điều chế PFT
            PftCodewordsBuffer(:,k)=wc_opt.'; % gán vectơ trọng số pha tối ưu wc_opt vào biến PftCodewordsBuffer ở cột thứ k
            PftGainBuffer(k)=array_gainpft; % gán giá trị độ lợi của tín hiệu điều chế PFT, dựa trên vectơ trọng số pha tối ưu, vào biến PftGainBuffer ở chỉ số k
        end

        %% tính toán các vectơ trọng số cho mỗi phương pháp (FC, NC và PFT) và sau đó nhân chúng với ma trận từ mã codebook tương ứng để tạo ra các vectơ trọng số cuối cùng cho mỗi phương pháp
        Product_mxg_DFT=prod(sqrt(FCGainBuffer)); % tích của căn bậc hai của các độ lợi tín hiệu từ mã codebook cho phương pháp FC
        MultiBeamFC_Orig=FCCodewordsBuffer*((Product_mxg_DFT./sqrt(FCGainBuffer))); % kết quả cuối cùng của vectơ trọng số cho phương pháp FC
        Product_mxg_NC=prod(sqrt(NCGainBuffer)); % tích của căn bậc hai của các độ lợi tín hiệu từ mã codebook cho phương pháp NC
        MultiBeamNC_Orig=NCCodewordsBuffer*((Product_mxg_NC./sqrt(NCGainBuffer))); % kết quả cuối cùng của vectơ trọng số cho phương pháp NC
        Product_mxg_Pft=prod(sqrt(PftGainBuffer)); % tích của căn bậc hai của các độ lợi tín hiệu từ mã codebook cho phương pháp PFT
        MultiBeamPft_Orig=PftCodewordsBuffer*((Product_mxg_Pft./sqrt(PftGainBuffer))); % kết quả cuối cùng của vectơ trọng số cho phương pháp PFT

        %% Điều chỉnh vectơ trọng số cuối cùng cho phương pháp FC để đảm bảo rằng không có trọng số nào bằng 0
        %MultiBeamFC_Orig=sum(FCCodewordsBuffer,2); % tính giá trị trung bình cộng của các vectơ trong FCCodewordsBuffer theo chiều cột, lưu kết quả vào MultiBeamFC_Orig
        record_zeroFC=find(MultiBeamFC_Orig==0); % tìm các vị trí trong MultiBeamFC_Orig mà có giá trị bằng 0 và lưu các chỉ mục đó vào mảng record_zeroFC
        MultiBeamFC_Orig(record_zeroFC)=exp(1j*2*pi*rand)/sqrt(N); %  thay thế các giá trị trong MultiBeamFC_Orig tại các vị trí có giá trị bằng 0 bằng một giá trị phức ngẫu nhiên được tạo ra từ hàm rand, sau đó chia cho căn bậc hai của N
        MultiBeamFCRIS=MultiBeamFC_Orig./abs(MultiBeamFC_Orig)/sqrt(N); % tạo ra một vectơ mới MultiBeamFCRIS bằng cách chia MultiBeamFC_Orig cho giá trị tuyệt đối của nó, sau đó chia kết quả cho căn bậc hai của N. Quá trình này làm cho vectơ MultiBeamFCRIS có các phần tử có độ lớn bằng 1 và được chuẩn hóa theo căn bậc hai của N
        MultiBeamFCAP=MultiBeamFC_Orig./max(abs(MultiBeamFC_Orig))/sqrt(N); % tạo ra một vectơ mới MultiBeamFCAP bằng cách chia MultiBeamFC_Orig cho giá trị tuyệt đối lớn nhất của nó, sau đó chia kết quả cho căn bậc hai của N. Quá trình này chuẩn hóa các phần tử của MultiBeamFC_Orig sao cho phần tử lớn nhất có độ lớn bằng 1, và sau đó chuẩn hóa kết quả theo căn bậc hai của N

        %% Chồng chập tia cận trường
        %MultiBeamNC_Orig=sum(NCCodewordsBuffer,2); % tính giá trị trung bình cộng của các cột trong ma trận NCCodewordsBuffer và lưu kết quả vào vectơ MultiBeamNC_Orig.
        record_zeroNC=find(MultiBeamNC_Orig==0); % tìm các vị trí trong vectơ MultiBeamNC_Orig có giá trị bằng 0 và lưu các vị trí đó vào vectơ record_zeroNC.
        MultiBeamNC_Orig(record_zeroNC)=exp(1j*2*pi*rand)/sqrt(N); % thay thế các giá trị 0 trong vectơ MultiBeamNC_Orig bằng một số phức ngẫu nhiên có độ lớn là 1 và góc pha ngẫu nhiên được chọn từ phạm vi từ 0 đến 2π để tránh chia cho 0 trong quá trình chuẩn hóa và tính toán về sau.
        MultiBeamNCRIS=MultiBeamNC_Orig./abs(MultiBeamNC_Orig)/sqrt(N); % chuẩn hóa vectơ MultiBeamNC_Orig bằng cách chia cho độ lớn tuyệt đối của nó, sau đó chia cho căn bậc hai của độ lớn N
        MultiBeamNCDig=MultiBeamNC_Orig/norm(MultiBeamNC_Orig); % tạo ra vectơ MultiBeamNCDig bằng cách chuẩn hóa MultiBeamNC_Orig theo chuẩn 2, có nghĩa là chia MultiBeamNC_Orig cho độ dài của nó (norm) tạo ra một vectơ mới có cùng hướng nhưng có độ lớn bằng 1.
        MultiBeamNCAP=MultiBeamNC_Orig./max(abs(MultiBeamNC_Orig))/sqrt(N);%Có thể thay đổi hệ số Biên độ và Pha. (Biên độ chạy trong khoảng [0,1])
        %
        %% Chồng chập tia Pft-BF ( kết hợp Precoding và Combining vào tạo tia )
        %MultiBeamPft_Orig=sum(PftCodewordsBuffer,2); %  tính giá trị trung bình cộng các cột của ma trận PftCodewordsBuffer, tạo ra vectơ MultiBeamPft_Orig.
        record_zeroPft=find(MultiBeamPft_Orig==0); % tìm các vị trí trong MultiBeamPft_Orig có giá trị bằng 0 và lưu các vị trí đó vào mảng record_zeroPft.
        MultiBeamPft_Orig(record_zeroPft)=exp(1j*2*pi*rand)/sqrt(N); % thay thế các giá trị trong MultiBeamPft_Orig tại các vị trí có giá trị bằng 0 bằng một giá trị ngẫu nhiên phức tạp được tính dựa trên hàm rand (chọn một số ngẫu nhiên từ phân phối đều trong khoảng từ 0 đến 1), sau đó chia cho căn bậc hai của N
        MultiBeamPftRIS=MultiBeamPft_Orig./abs(MultiBeamPft_Orig)/sqrt(N); % tính toán MultiBeamPftRIS bằng cách chia MultiBeamPft_Orig cho giá trị tuyệt đối của nó, sau đó chia cho căn bậc hai của N, chuẩn hóa vector MultiBeamPft_Orig và đưa về cùng một phạm vi với các giá trị nằm trong khoảng từ 0 đến 1.
        MultiBeamPftDig=MultiBeamPft_Orig/norm(MultiBeamPft_Orig); % tạo MultiBeamPftDig bằng cách chuẩn hóa MultiBeamPft_Orig bằng cách chia nó cho norm (chuẩn) của chính nó. Điều này sẽ biến MultiBeamPftDig thành một vector có độ dài (độ lớn) bằng 1, giữ nguyên hướng của vector gốc.
        MultiBeamRand=exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N); % tạo một vector MultiBeamRand gồm các phần tử phức hợp, mỗi phần tử có giá trị là một số phức hợp được tạo ra ngẫu nhiên trong phạm vi [-1, 1] (tương đương với góc từ -180 đến 180 độ). Sau đó, nó được chuẩn hóa bằng cách chia cho căn bậc hai của độ dài của vector để đảm bảo rằng norm của vector là 1.

        %
        %% Tính độ lợi đa tia
        G_MultiBeam_FFSuperpose=abs(MultiBeamFCRIS.'*GG).^2;%size:1 * K (RIS BF) % tính toán giá trị của mảng G_MultiBeam_FFSuperpose. Đầu tiên, nó nhân tích vô hướng của vector cột MultiBeamFCRIS và ma trận GG, sau đó lấy giá trị tuyệt đối bình phương của kết quả. 
        G_MultiBeam_NFSuperpose=abs(MultiBeamNCRIS.'*GG).^2;%(RIS BF) % tính toán giá trị của mảng G_MultiBeam_NFSuperpose. Tương tự như trước, nó nhân tích vô hướng của vector cột MultiBeamNCRIS và ma trận GG, sau đó lấy giá trị tuyệt đối bình phương của kết quả.
        G_MultiBeam_PftSuperpose=abs(MultiBeamPftDig.'*GG).^2;%(Dig BF) % tính toán giá trị của mảng G_MultiBeam_PftSuperpose. Nó thực hiện phép nhân tích vô hướng của vector hàng MultiBeamPftDig và ma trận GG, sau đó lấy giá trị tuyệt đối bình phương của kết quả. 
        G_MultiBeam_FFSuperpose_AP=abs(MultiBeamFCAP.'*GG).^2; % tính toán giá trị của mảng G_MultiBeam_FFSuperpose_AP. Nó thực hiện phép nhân tích vô hướng của vector hàng MultiBeamFCAP và ma trận GG, sau đó lấy giá trị tuyệt đối bình phương của kết quả.
        G_MultiBeam_NFSuperpose_AP=abs(MultiBeamNCAP.'*GG).^2; % tính toán giá trị của mảng G_MultiBeam_NFSuperpose_AP. Nó thực hiện phép nhân tích vô hướng của vector hàng MultiBeamNCAP và ma trận GG, sau đó lấy giá trị tuyệt đối bình phương của kết quả.
        G_MultiBeam_Rand=abs(MultiBeamRand.'*GG).^2; % tính toán giá trị của mảng G_MultiBeam_Rand. Nó thực hiện phép nhân tích vô hướng của vector hàng MultiBeamRand và ma trận GG, sau đó lấy giá trị tuyệt đối bình phương của kết quả.

        %% Thuật toán MM
        loss3 = []; % tạo mảng lưu giá trị thất thoát (loss)
        eta=1; % hệ số điều chỉnh vector trọng số v_abs3
        v_abs3=max(sqrt(G_MultiBeam_PftSuperpose))*ones(1,num_K)'*eta; % vector cột có kishc thước num_K x 1 với mỗi phần tử được đặt bằng giá trị lớn nhất trong dãy các giá trị bình phương của G_MultiBeam_PftSuperpose nhân với eta.
        w3 = exp(1j*1*pi*(2*rand(N,1)-1))/sqrt(N); % vector hàng với các phần tử được khởi tạo ngẫu nhiên theo phân phối đều và được chuẩn hóa.
        
        v_phase=exp(1j*1*pi*(2*rand(num_K,1)-1)); % vector cột có kích thước num_K x 1, với các phần tử được khởi tạo ngẫu nhiên theo phân phối đều.
        v3 = v_abs3.*v_phase; % vector cột có kích thước num_K x 1, kết hợp giữa v_abs3 và v_phase.
        
        A=conj(NCCodewordsBuffer).';%dùng mã cận trường k*n
        A=diag(sqrt(NCGainBuffer))*conj(NCCodewordsBuffer).'; % ma trận có kích thước k x n (với k là số lượng người dùng và n là số lượng anten), được tạo ra từ cơ sở anten của người dùng (NCCodewordsBuffer) và điều chỉnh bởi các hệ số truyền (NCGainBuffer) cung cấp thông tin về cơ sở anten của người dùng trong mô hình kênh.
        %A=GG.';
        lambda = max(eig(A'*A));
        max_iter=200; % số lần lặp tối đa

        % Vòng lặp
        for i=1:max_iter
            f = A*w3; % cập nhật hệ số pha v_phase
            v_phase = exp(1j * angle(f)); % cập nhật các hệ số pha
            v3 = v_abs3 .* v_phase; % cập nhật hệ số w3
            
            w3_0iter=w3; % sao chép giá trị của vectơ w3 vào w3_0iter để sử dụng làm giá trị khởi tạo cho việc cập nhật vectơ w3 trong quá trình lặp.
            for j_MM=1:10
                temp = A' * v3 - A' * A * w3_0iter + lambda * w3_0iter ;
                w3_0iter = exp( 1j*angle( temp ))/sqrt(N);
                
            end
            w3=w3_0iter;
            loss3 = [loss3, norm(v3 - A * w3, 2)]; % tính toán và cập nhật giá trị của loss3, là giá trị trung bình cộng các norm 2 của sự sai khác giữa v3 và A * w3, và nó được sử dụng để theo dõi sự hội tụ của quá trình cập nhật.
        end

        %% tính độ lợi tia theo thuật toán MM
        G_MultiBeam_MM=abs(w3.'*GG).^2; %  tính toán giá trị của G_MultiBeam_MM, là giá trị bình phương của norm 2 của tích vô hướng giữa w3 và GG.
        G_MM=abs(w3.'*GG).^2; % tính toán giá trị của G_MM, là giá trị bình phương của norm 2 của tích vô hướng giữa w3 và GG.

        %% tính tốc độ truyền cho từng người dùng
        R_FF_RIS=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose.');%Size K*1 % hiệu quả phổ của hệ thống khi sử dụng phương pháp Beamforming tại RIS (RIS BF).
        R_NF_RIS=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose.'); % hiệu quả phổ của hệ thống khi sử dụng phương pháp Near-field Beamforming tại RIS (NF RIS).
        R_FF_AP=log2(1+SNR_linear.*G_MultiBeam_FFSuperpose_AP.'); % Sự hiệu quả phổ của hệ thống khi sử dụng Beamforming tại AP (AP BF).
        R_NF_AP=log2(1+SNR_linear.*G_MultiBeam_NFSuperpose_AP.'); % Sự hiệu quả phổ của hệ thống khi sử dụng Near-field Beamforming tại AP (NF AP).
        R_MM=log2(1+SNR_linear.*G_MultiBeam_MM.'); % Sự hiệu quả phổ của hệ thống khi sử dụng phương pháp Multiple Measurement Vectors (MM).
        R_Rand=log2(1+SNR_linear.*G_MultiBeam_Rand.'); % Sự hiệu quả phổ của hệ thống khi sử dụng phương pháp Random (AP BF).
         
        %% Calculate Sum-Rate & MinRate
        SumR_FF_RIS=sum(R_FF_RIS); % Tổng tỷ lệ thông tin cho RIS Beamforming.
        SumR_NF_RIS=sum(R_NF_RIS); % Tổng tỷ lệ thông tin cho RIS Near-field Beamforming.
        SumR_FF_AP=sum(R_FF_AP); % Tổng tỷ lệ thông tin cho AP Beamforming.
        SumR_NF_AP=sum(R_NF_AP); % Tổng tỷ lệ thông tin cho AP Near-field Beamforming.
        SumR_MM=sum(R_MM); % Tổng tỷ lệ thông tin cho Multiple Measurement Vectors.
        SumR_Rand=sum(R_Rand); % Tổng tỷ lệ thông tin cho kỹ thuật ngẫu nhiên.
        
        MinR_FF_RIS=min(R_FF_RIS); % Giá trị tối thiểu của tỷ lệ thông tin cho RIS Beamforming.
        MinR_NF_RIS=min(R_NF_RIS); % Giá trị tối thiểu tỷ lệ thông tin cho RIS Near-field Beamforming.
        MinR_FF_AP=min(R_FF_AP); % Giá trị tối thiểu tỷ lệ thông tin cho AP Beamforming.
        MinR_NF_AP=min(R_NF_AP); % Giá trị tối thiểu tỷ lệ thông tin cho AP Near-field Beamforming.
        MinR_MM=min(R_MM); % Giá trị tối thiểu tỷ lệ thông tin cho Multiple Measurement Vectors.
        MinR_Rand=min(R_Rand); % Giá trị tối thiểu tỷ lệ thông tin cho kỹ thuật ngẫu nhiên.

        % Tích hợp giá trị trung bình của tỷ lệ thông tin và giá trị tối thiểu của nó qua mỗi vòng lặp của quá trình mô phỏng hoặc tính toán. 
        % Giúp ổn định kết quả và cung cấp ước lượng trung bình của hiệu suất của từng phương pháp kết hợp sóng.
        % Đặc biệt, việc chia cho số lượng lặp (ITER) giúp xác định giá trị trung bình dựa trên nhiều lần thử nghiệm.
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
    %% Lưu thông số tỷ lệ sau khi kết thúc vòng lặp
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
    
            
        Record_SumR_Rand(idx_Bigfor)=temp_SumR_Rand;
        
        
        Record_MinR_Rand(idx_Bigfor)=temp_MinR_Rand;
    
end

xx=Bigfor_list;

%% Lưu vào file
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

%% Vẽ bằng số liệu có sẵn.
%Muốn lấy Fig liền thì xóa comment rồi chạy từ đây đến hết code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('xx.mat');
% load('Record_MinR_FF_RIS.mat');
% load('Record_MinR_NF_RIS.mat');
% load('Record_MinR_FF_AP.mat');
% load('Record_MinR_NF_AP.mat');
% load('Record_MinR_MM.mat');
% load('Record_MinR_Rand.mat');
% %-------------------------------%
% load('Record_SumR_FF_RIS.mat');
% load('Record_SumR_NF_RIS.mat');
% load('Record_SumR_FF_AP.mat');
% load('Record_SumR_NF_AP.mat');
% load('Record_SumR_MM.mat');
% load('Record_SumR_Rand.mat');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% END-Nếu muốn lấy Fig liền với số liệu sẵn có.


label_fontsize=13;
legend_fontsize=11;
linewidth=1.5;

%% Vẽ đồ thị tỷ số truyền tối thiểu
figure; hold on; box on; grid on;
p1 = plot(xx,Record_MinR_MM, 'Color', [246,83,20]./255,  'Linestyle', '-',  'Marker', '^'   , 'Linewidth', linewidth);
p2 = plot(xx,Record_MinR_FF_RIS, 'Color', [124,187,0]./255,  'Linestyle', '-',  'Marker', 's' , 'Linewidth', linewidth);
p3 = plot(xx,Record_MinR_NF_RIS, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', 'd'  , 'Linewidth', linewidth);
% p4 = plot(xx,Record_MinR_FF_AP, 'Color', [124,187,0]./255, 'Linestyle', '-',  'Marker', '<'    , 'Linewidth', linewidth);
% p5 = plot(xx,Record_MinR_NF_AP, 'Color', [0,161,241]./255,  'Linestyle', '-',  'Marker', '>'   , 'Linewidth', linewidth);
p6 = plot(xx,Record_MinR_Rand, 'Color', 255.*[0.49,0.18,0.56]./255,  'Linestyle', '--',  'Marker', 'x'  , 'Linewidth', linewidth);

%set(gca,'XDir','reverse')
legend([p1,p3,p2,p6], {'Proposed beam design', 'NFCB-based beam design','FFCB-based beam design','Random phase design'},'fontsize',legend_fontsize)
xlabel('SNR [dB]','fontsize',label_fontsize)
ylabel('Minimal rate [bits/s/Hz]','fontsize',label_fontsize)

toc
