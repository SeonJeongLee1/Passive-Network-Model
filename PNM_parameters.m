clear
V = 3.5; % set voltage
V2 = -1.0; % reset voltag
step = 0.01; % step voltage
alpha=  1; 
beta = 0;
m = 10; % Unified node 'm' for all layers
n = 40;
Vset = 1.75/m; % Setting Vset and Vreset voltages across each partial resistance
Vreset = 0.82/m;
min_black = 2;

% Adjustable probability value
fail1 = 0; % Probability of failure during the reset process
fail2 = 0; % Probability of failure during the set process

aaa = 0; % Standard deviation
bbb = 20; % Standard deviation of Vset voltage
ccc = 20; % Standard deviation of Vreset voltage

% Defect probability for Top and Bottom 
percent_defect = 20; %Top defect [%]
percent_defect3 = 1; %Bottom defect [%]
num_splits=m-1; % To achieve a gradual distribution of LRS (defect) values across rows, the 'm' node is divided

% Distribution of Vset and Vreset voltages
Vset_change = Vset/(100/bbb);
Vreset_change = Vreset/(100/ccc);
Vset=Vset+Vset_change*randn(1,1);
Vreset=Vreset+Vreset_change*randn(1,1);
[Vset_h, Vreset_h] = variation_Voltage(Vset, Vreset, Vset_change, Vreset_change, (n-1)*(m-2));
[Vset_v, Vreset_v] = variation_Voltage(Vset, Vreset, Vset_change, Vreset_change, n*(m-1));
[Vset_d1, Vreset_d1] = variation_Voltage(Vset, Vreset, Vset_change, Vreset_change, (n-1)*(m-1));
[Vset_d2, Vreset_d2] = variation_Voltage(Vset, Vreset, Vset_change, Vreset_change, (n-1)*(m-1));

% Defining resistance values
a1=7.07*10^-6; % Resistance associated with Schottky emission
a2=6.33;
a3=-6.47;
HRS = 2309538.186*n/(m*0.414); % Resistance associated with Ohmic current
LRS = 239*n/(m*0.414);

% High resistance function for Schottky emission
R_h=@(VV) VV/(a1*exp(a2*sqrt(VV)+a3)); % HRS
Rn_h=@(VV) R_h(VV)*n/(m*0.414);
Rh_h_A = @(VV) Rn_h(VV);  % Average high resistance (Rh) 
Rv_h_A = Rn_h;  % Average high resistance (Rv) 
Rd_h_A = @(VV) Rh_h_A(VV)*sqrt(2);  % Average high resistance (Rd) 
layer_Rh_s_A=@(VV) Rh_h_A(VV)/(100/aaa);    % Adjustment of standard deviation of high resistance (Rh) 
layer_Rv_s_A=@(VV) Rv_h_A(VV)/(100/aaa);   % Adjustment of standard deviation of high resistance (Rv)
layer_Rd_s_A=@(VV) Rd_h_A(VV)/(100/aaa);   % Adjustment of standard deviation of high resistance (Rd)

% Setting individual resistance values using the Schottky emission model
Rh_h_B = HRS;  % Average high resistance (Rh)
Rv_h_B = HRS;  % Average high resistance (Rv)
Rd_h_B = Rh_h_B*sqrt(2);  % Average high resistance (Rd)
layer_Rh_s_B=Rh_h_B/(100/aaa);   % Adjustment of standard deviation for high resistance Rh
layer_Rv_s_B=Rv_h_B/(100/aaa);   % Adjustment of standard deviation for high resistance Rv
layer_Rd_s_B=Rd_h_B/(100/aaa);   % Adjustment of standard deviation for high resistance Rd

% Low resistance function for ohmic behavior
R_l=@(VV) VV/(a1*exp(a2*sqrt(VV)+a3)); % LRS
Rn_l=@(VV) R_l(VV)*n/(m*0.414);
Rh_l_A = @(VV) Rn_l(VV);  % Average of low resistance Rh
Rv_l_A = Rn_l;  % Average of low resistance Rv
Rd_l_A = @(VV) Rh_l_A(VV)*sqrt(2);  % Average of low resistance Rd
layer_Rh_ss_A=@(VV) Rh_l_A(VV)/(100/aaa);   % Adjustment of standard deviation for low resistance Rh
layer_Rv_ss_A=@(VV) Rv_l_A(VV)/(100/aaa);   % Adjustment of standard deviation for low resistance Rv
layer_Rd_ss_A=@(VV) Rd_l_A(VV)/(100/aaa);   % Adjustment of standard deviation for low resistance Rd

% Setting individual resistance values using ohmic behavior
Rh_l_B = LRS;  % Average of low resistance Rh
Rv_l_B = LRS;  % Average of low resistance Rv
Rd_l_B = Rh_l_B*sqrt(2);  % Average of low resistance Rd
layer_Rh_ss_B=Rh_l_B/(100/aaa);   % Adjustment of standard deviation for low resistance Rh
layer_Rv_ss_B=Rv_l_B/(100/aaa);   % Adjustment of standard deviation for low resistance Rv
layer_Rd_ss_B=Rd_l_B/(100/aaa);   % Adjustment of standard deviation for low resistance Rd

voltage1 = 0: step :V;
voltage2 = V : -step :0;
voltage3 = 0: -step : V2;
voltage4 = V2: step : 0;

% Speed improvement through preallocation
Rt1=zeros(1,length(voltage1)); Rt2=zeros(1,length(voltage2)); Rt3=zeros(1,length(voltage3)); Rt4=zeros(1,length(voltage4));
layer_Rh_hh_A(1:(n-1)*(m-2))={0}; 
layer_Rv_hh_A(1:n*(m-1))={0}; 
layer_Rd_1hh_A(1:(n-1)*(m-1))={0}; 
layer_Rd_2hh_A(1:(n-1)*(m-1))={0}; 

layer_Rh_hh_B=zeros(1,(n-1)*(m-2)); 
layer_Rv_hh_B=zeros(1,n*(m-1)); 
layer_Rd_1hh_B=zeros(1,(n-1)*(m-1)); 
layer_Rd_2hh_B=zeros(1,(n-1)*(m-1)); 

layer_Rh_ll_A(1:(n-1)*(m-2))={0}; 
layer_Rv_ll_A(1:n*(m-1))={0}; 
layer_Rd_1ll_A(1:(n-1)*(m-1))={0}; 
layer_Rd_2ll_A(1:(n-1)*(m-1))={0}; 

layer_Rh_ll_B=zeros(1,(n-1)*(m-2)); 
layer_Rv_ll_B=zeros(1,n*(m-1)); 
layer_Rd_1ll_B=zeros(1,(n-1)*(m-1)); 
layer_Rd_2ll_B=zeros(1,(n-1)*(m-1)); 

%-----------------------------------------------schottky emission_HRS
%%%% Set horizontal resistance values
layer_Rh_hh_A_ran = randn(1,(n-1)*(m-2));
for i = 1 : (n-1)*(m-2)                          % Set up the Rh function of high resistance reflecting the standard deviation
    layer_Rh_hh_A{i} = @(VV) layer_Rh_s_A(VV)*layer_Rh_hh_A_ran(i)+Rh_h_A(VV);
end

Rh_hh_A = layer_Rh_hh_A;

%%%%% Set vertical resistance values
layer_Rv_hh_A_ran = randn(1,n*(m-1));
for i = 1 : n*(m-1)                         % Set up the Rv function for high resistance reflecting the standard deviation
    layer_Rv_hh_A{i} = @(VV) layer_Rv_s_A(VV)*layer_Rv_hh_A_ran(i)+Rv_h_A(VV);
end

Rv_hh_A = layer_Rv_hh_A;
%%%%%

%%%%% Set diagonal resistance values
layer_Rd_1hh_A_ran = randn(1,(n-1)*(m-1));
layer_Rd_2hh_A_ran = randn(1,(n-1)*(m-1));
for i = 1 : (n-1)*(m-1)                     % Set up the Rd function for high resistance reflecting the standard deviation
    layer_Rd_1hh_A{i} = @(VV) layer_Rd_s_A(VV)*layer_Rd_1hh_A_ran(i)+Rd_h_A(VV);
    layer_Rd_2hh_A{i} = @(VV) layer_Rd_s_A(VV)*layer_Rd_2hh_A_ran(i)+Rd_h_A(VV);
end

Rd_1hh_A = layer_Rd_1hh_A;
Rd_2hh_A = layer_Rd_2hh_A;

%-----------------------------------------------ohmic current_HRS
layer_Rh_hh_B_ran=randn(1,(n-1)*(m-2));
%%%% Set horizontal resistance values
for i = 1 : (n-1)*(m-2)                         % Set up the Rh function of high resistance reflecting the standard deviation
   layer_Rh_hh_B(i) = layer_Rh_s_B*layer_Rh_hh_B_ran(i)+Rh_h_B;
end

Rh_hh_B = layer_Rh_hh_B;

layer_Rv_hh_B_ran=randn(1,n*(m-1));
for i = 1 : n*(m-1)                         % Set up the Rv function for high resistance reflecting the standard deviation
   layer_Rv_hh_B(i) = layer_Rv_s_B*layer_Rv_hh_B_ran(i)+Rv_h_B;
end

Rv_hh_B = layer_Rv_hh_B;
%%%%%

%%%%% Set diagonal resistance values
layer_Rd_1hh_B_ran=randn(1,(n-1)*(m-1));
layer_Rd_2hh_B_ran=randn(1,(n-1)*(m-1));
for i = 1 : (n-1)*(m-1)                     % Set up the Rd function for high resistance reflecting the standard deviation
   layer_Rd_1hh_B(i) = layer_Rd_s_B*layer_Rd_1hh_B_ran(i)+Rd_h_B;
   layer_Rd_2hh_B(i) = layer_Rd_s_B*layer_Rd_2hh_B_ran(i)+Rd_h_B;
end

Rd_1hh_B = layer_Rd_1hh_B;
Rd_2hh_B = layer_Rd_2hh_B;

%%%%%LRS
%-----------------------------------------------schottky emission_LRS
%%%% Set horizontal resistance values
layer_Rh_ll_A_ran = randn(1,(n-1)*(m-2));
for i = 1 : (n-1)*(m-2)                         % Set up the Rh function of low resistance reflecting the standard deviation
    layer_Rh_ll_A{i} = @(VV) layer_Rh_ss_A(VV)*layer_Rh_ll_A_ran(i)+Rh_l_A(VV);
end

Rh_ll_A = layer_Rh_ll_A;

%%%%% Set vertical resistance values
layer_Rv_ll_A_ran = randn(1,n*(m-1));
for i = 1 : n*(m-1)                         % Set up the Rv function for low resistance reflecting the standard deviation
    layer_Rv_ll_A{i} = @(VV) layer_Rv_ss_A(VV)*layer_Rv_ll_A_ran(i)+Rv_l_A(VV);
end

Rv_ll_A = layer_Rv_ll_A;
%%%%%

%%%%% Set diagonal resistance values
layer_Rd_1ll_A_ran = randn(1,(n-1)*(m-1));
layer_Rd_2ll_A_ran = randn(1,(n-1)*(m-1));
for i = 1 : (n-1)*(m-1)                     % Set up the Rd function for low resistance reflecting the standard deviation
    layer_Rd_1ll_A{i} = @(VV) layer_Rd_ss_A(VV)*layer_Rd_1ll_A_ran(i)+Rd_l_A(VV);
    layer_Rd_2ll_A{i} = @(VV) layer_Rd_ss_A(VV)*layer_Rd_2ll_A_ran(i)+Rd_h_A(VV);
end

Rd_1ll_A = layer_Rd_1ll_A;
Rd_2ll_A = layer_Rd_2ll_A;

%-----------------------------------------------ohmic current_LRS
%%%% Set horizontal resistance values
layer_Rh_ll_B_ran =randn(1,(n-1)*(m-2));
for i = 1 : (n-1)*(m-2)                         % Set up the Rh function of low resistance reflecting the standard deviation
    layer_Rh_ll_B(i) = layer_Rh_ss_B*layer_Rh_ll_B_ran(i)+Rh_l_B;
end

Rh_ll_B = layer_Rh_ll_B;

%%%%% Set vertical resistance values
layer_Rv_ll_B_ran =randn(1,n*(m-1));
for i = 1 : n*(m-1)                         % Set up the Rv function for low resistance reflecting the standard deviation
    layer_Rv_ll_B(i) = layer_Rv_ss_B*layer_Rv_ll_B_ran(i)+Rv_l_B;
end

Rv_ll_B = layer_Rv_ll_B;
%%%%%

%%%%% Set diagonal resistance values
layer_Rd_1ll_B_ran=randn(1,(n-1)*(m-1));
layer_Rd_2ll_B_ran = randn(1,(n-1)*(m-1));
for i = 1 : (n-1)*(m-1)                     % Set up the Rd function for low resistance reflecting the standard deviation
    layer_Rd_1ll_B(i) = layer_Rd_ss_B*layer_Rd_1ll_B_ran(i)+Rd_l_B;
    layer_Rd_2ll_B(i) = layer_Rd_ss_B*layer_Rd_2ll_B_ran(i)+Rd_l_B;
end

Rd_1ll_B = layer_Rd_1ll_B;
Rd_2ll_B = layer_Rd_2ll_B;

%%%%% Save the initial resistance state
Rh_cell_A=Rv_hh_A; Rv_cell_A=Rv_hh_A; 
Rd_1cell_A=Rd_1hh_A; Rd_2cell_A=Rd_2hh_A;

Rh_cell2_A=Rv_hh_A; Rv_cell2_A=Rv_hh_A; 
Rd_1cell2_A=Rd_1hh_A; Rd_2cell2_A=Rd_2hh_A;

Rh_cell_Al=Rv_ll_A; Rv_cell_Al=Rv_ll_A; 
Rd_1cell_Al=Rd_1ll_A; Rd_2cell_Al=Rd_1ll_A;

Rh_cell2_Al=Rv_ll_A; Rv_cell2_Al=Rv_ll_A; 
Rd_1cell2_Al=Rd_1ll_A; Rd_2cell2_Al=Rd_1ll_A;

%%%%% Set the initial voltage input for the Schottky emission function 
for i=1:(n-1)*(m-2)
    Rh_cell_A{i}=Rh_cell2_A{i}(voltage1(1,2));
    Rh_A(i)=Rh_cell_A{i};
end

for i=1:n*(m-1)
    Rv_cell_A{i}=Rv_cell2_A{i}(voltage1(1,2));
    Rv_A(i)=Rv_cell_A{i};
end

for i=1:(n-1)*(m-1)
    Rd_1cell_A{i}=Rd_1cell2_A{i}(voltage1(1,2));
    Rd_1_A(i)=Rd_1cell_A{i};
end

for i=1:(n-1)*(m-1)
    Rd_2cell_A{i}=Rd_2cell2_A{i}(voltage1(1,2));
    Rd_2_A(i)=Rd_2cell_A{i};
end
%%%%% Set the initial voltage input for the Ohmic conduction function 
for i=1:(n-1)*(m-2)
    Rh_cell_Al{i}=Rh_cell2_Al{i}(voltage1(1,2));
    Rh_Al(i)=Rh_cell_Al{i};
end

for i=1:n*(m-1)
    Rv_cell_Al{i}=Rv_cell2_Al{i}(voltage1(1,2));
    Rv_Al(i)=Rv_cell_Al{i};
end

for i=1:(n-1)*(m-1)
    Rd_1cell_Al{i}=Rd_1cell2_Al{i}(voltage1(1,2));
    Rd_1_Al(i)=Rd_1cell_Al{i};
end

for i=1:(n-1)*(m-1)
    Rd_2cell_Al{i}=Rd_2cell2_Al{i}(voltage1(1,2));
    Rd_2_Al(i)=Rd_2cell_Al{i};
end

Rh_hh=alpha*Rh_A+(1-alpha)*Rh_h_B;
Rv_hh=alpha*Rv_A+(1-alpha)*Rv_h_B;
Rd_1hh=alpha*Rd_1_A+(1-alpha)*Rd_h_B;
Rd_2hh=alpha*Rd_2_A+(1-alpha)*Rd_h_B;

pre_Rh_hh=Rh_hh; pre_Rv_hh=Rv_hh; pre_Rd_1hh=Rd_1hh; pre_Rd_2hh=Rd_2hh;

Rh_ll=beta*Rh_Al+(1-beta)*Rh_l_B;
Rv_ll=beta*Rv_Al+(1-beta)*Rv_l_B;
Rd_1ll=beta*Rd_1_Al+(1-beta)*Rd_l_B;
Rd_2ll=beta*Rd_2_Al+(1-beta)*Rd_l_B;

%%%%% Set up the code to distribute the defects
increment_Rh = n-1;
increment_Rv = n;
increment_Rd = n-1;
Rh_ranges = [(1:increment_Rh:length(Rh_hh))', min((increment_Rh+1:increment_Rh:length(Rh_hh) + increment_Rh)', length(Rh_hh))];
Rv_ranges = [(1:increment_Rv:length(Rv_hh))', min((increment_Rv+1:increment_Rv:length(Rv_hh) + increment_Rv)', length(Rv_hh))];
Rd_1ranges = [(1:increment_Rd:length(Rd_1hh))', min((increment_Rd+1:increment_Rd:length(Rd_1hh) + increment_Rd)', length(Rd_1hh))];
Rd_2ranges = [(1:increment_Rd:length(Rd_2hh))', min((increment_Rd+1:increment_Rd:length(Rd_2hh) + increment_Rd)', length(Rd_2hh))];

%%%% Create a probability vector
step2 = (percent_defect - percent_defect3) / (num_splits - 1);
prob_vec = percent_defect:-step2:percent_defect3;

for k = 1:m-1
    % Combine the entire index range
    if k == 1
        all_ranges = [Rv_ranges(k, :), Rd_1ranges(k, :), Rd_2ranges(k, :)];
    else
        all_ranges = [Rh_ranges(k-1, :), Rv_ranges(k, :), Rd_1ranges(k, :), Rd_2ranges(k, :)];
    end

    % Calculate the total number of elements
    total_elements = sum(all_ranges(2:2:end) - all_ranges(1:2:end));

    % Calculate the number of elements based on probability for the entire set of elements
    num_defects = round(total_elements * prob_vec(k) / 100);

    % Shuffle all the indexes
    all_idx = randperm(total_elements);

    % Select the 'num_defects' number of indexes
    selected_idx = all_idx(1:num_defects);

    % Convert the selected indexes to the indexes of each vector and change the values
    for i = 1:length(selected_idx)
        idx = selected_idx(i);

        % Initialize lengths
        if k > 1
            Rh_length = Rh_ranges(k-1, 2) - Rh_ranges(k-1, 1);
        else
            Rh_length = 0;
        end
        Rv_length = Rv_ranges(k, 2) - Rv_ranges(k, 1);
        Rd_1_length = Rd_1ranges(k, 2) - Rd_1ranges(k, 1);

        if k > 1 && idx <= Rh_length
            Rh_hh(Rh_ranges(k-1, 1) + idx - 1) = Rh_ll(Rh_ranges(k-1, 1) + idx - 1);
        elseif idx <= Rh_length + Rv_length
            idx = idx - Rh_length;
            Rv_hh(Rv_ranges(k, 1) + idx - 1) = Rv_ll(Rv_ranges(k, 1) + idx - 1);
        elseif idx <= Rh_length + Rv_length + Rd_1_length 
            idx = idx - Rh_length - Rv_length;
            Rd_1hh(Rd_1ranges(k, 1) + idx - 1) = Rd_1ll(Rd_1ranges(k, 1) + idx - 1);
        else 
            idx = idx - Rh_length - Rv_length - Rd_1_length;
            Rd_2hh(Rd_2ranges(k, 1) + idx - 1) = Rd_2ll(Rd_2ranges(k, 1) + idx - 1);
        end
    end
end

Rh=Rh_hh; Rv=Rv_hh; % Set the initial high resistance state
Rd_1=Rd_1hh;
Rd_2=Rd_2hh;

A_Rh=reshape(Rh, n-1, m-2)';
A_Rv=reshape(Rv, n,m-1)';
A_Rd1=reshape(Rd_1, n-1,m-1)';
A_Rd2=reshape(Rd_2, n-1,m-1)';

%% Create a resistor image
jj=1;
max=30000;
min=3000;

Rh=Rh'; Rv=Rv'; Rd_1=Rd_1'; Rd_2=Rd_2';
Rh_ll=Rh_ll'; Rv_ll=Rv_ll'; Rd_1ll=Rd_1ll'; Rd_2ll=Rd_2ll';

A_image_R=zeros(2*m-3,2*n-1);

%Rv value
k=1;
for i = 1 : m-1
    for j = 1 : n
        A_image_R(2*i-1,2*j-1) = Rv(k);
        k=k+1;
    end
end

%Rh value
k=1;
for i = 1 : m-2
    for j = 1 : n-1
        A_image_R(2*i,2*j) = Rh(k);
        k=k+1;
    end
end

vir_Rh=zeros(m-2,n-1); vir_Rh_ll=zeros(m-2,n-1); vir_Rv=zeros(m-1,n); vir_Rv_ll=zeros(m-1,n); vir_Rd_1=zeros(m-1,n-1); vir_Rd_1ll=zeros(m-1,n-1); vir_Rd_2=zeros(m-1,n-1); vir_Rd_2ll=zeros(m-1,n-1);

k=1;
for i = 1:m-2
    for j=1:n-1
        vir_Rh(i,j)=Rh(k,1);
        vir_Rh_ll(i,j)=Rh_ll(k,jj);
        k=k+1;
    end
end

k=1;
for i = 1:m-1
    for j=1:n
        vir_Rv(i,j)=Rv(k,1);
        vir_Rv_ll(i,j)=Rv_ll(k,jj);
        k=k+1;
    end
end

k=1;
for i = 1:m-1
    for j=1:n-1
        vir_Rd_1(i,j)=Rd_1(k,1);
        vir_Rd_1ll(i,j)=Rd_1ll(k,jj);
        k=k+1;
    end
end

k=1;
for i = 1:m-1
    for j=1:n-1
        vir_Rd_2(i,j)=Rd_2(k,1);
        vir_Rd_2ll(i,j)=Rd_2ll(k,jj);
        k=k+1;
    end
end

vir_Rh = [zeros(size(vir_Rh, 1), 1) vir_Rh zeros(size(vir_Rh, 1), 1)];
vir_Rh_ll = [ones(size(vir_Rh_ll, 1), 1) vir_Rh_ll ones(size(vir_Rh_ll, 1), 1)];
vir_Rd_1 = [zeros(size(vir_Rd_1, 1), 1) vir_Rd_1 zeros(size(vir_Rd_1, 1), 1)];
vir_Rd_1ll=[ones(size(vir_Rd_1ll, 1), 1) vir_Rd_1ll ones(size(vir_Rd_1ll, 1), 1)];
vir_Rd_2 = [zeros(size(vir_Rd_2, 1), 1) vir_Rd_2 zeros(size(vir_Rd_2, 1), 1)];
vir_Rd_2ll=[ones(size(vir_Rd_2ll, 1), 1) vir_Rd_2ll ones(size(vir_Rd_2ll, 1), 1)];

k=1;
for i = 1 : m-2
    for j = 1 : n
        A_image_R(2*i,2*j-1) = max;
        k=k+1;
    end
end

for i = 1 : m-2
    for j = 1 : n
        conditions = [vir_Rh(i,j)==vir_Rh_ll(i,j), vir_Rh(i,j+1)==vir_Rh_ll(i,j+1), vir_Rv(i,j)==vir_Rv_ll(i,j), vir_Rv(i+1,j)==vir_Rv_ll(i+1,j), vir_Rd_1(i,j)==vir_Rd_1ll(i,j), vir_Rd_1(i+1,j+1)==vir_Rd_1ll(i+1,j+1), vir_Rd_2(i,j+1)==vir_Rd_2ll(i,j+1), vir_Rd_2(i+1,j)==vir_Rd_2ll(i+1,j)];      
        % Check if at least two of the conditions are true
        if sum(conditions) >= min_black
            A_image_R(2*i,2*j-1) = min;
        end
    end
end

Rd_mean=zeros(size(Rd_1)); 
for i = 1 : size(Rd_1,1)
    for j =1:size(Rd_1,2)
        if Rd_1(i,j)<Rd_2(i,j)
            Rd_mean(i,j)=Rd_1(i,j);
        else
            Rd_mean(i,j)=Rd_2(i,j);
        end
    end
end

k=1;
for i = 1 : m-1
    for j = 1 : n-1
        A_image_R(2*i-1,2*j) = Rd_mean(k);
        k=k+1;
    end
end

clims= [min max];
imagesc(A_image_R, clims)
colorbar
colormap(gray)