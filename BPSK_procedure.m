clear all;
close all;

% �������� ������ BPSK ������������������
size = 100;
% ��������� ������� �������
fc = 10e6;
% ������� �������������
fs = 50e6;
% ����������� ������������
L = 50;  % �� ���� ��� �������� ��������� ����� �� 1 ������ ���������

SIGNAL_DURATION_micro = (size*L/fs)*1e6; % ������������ ������������� ��������� � �������������

%% ��������� �������


% ������������ ������� ������������������
RANDOMDATA = randi([0 1], 1, 10000);  % �������� ���������������� �������� ������
DATA = RANDOMDATA(1:size);            % ��������� ������� ���������� �����������

% ������������� BPSK ���������

for i = 1:length(DATA)
  if (DATA(i) == 1)
    MOD_DATA(i) = 1;
  else                    % ���� ���������� -1, � ������� �������� ���������
    MOD_DATA(i) = -1;
  end
end

% MOD_DATA_FFT = fft(MOD_DATA);
% 
% figure;
% subplot(2, 1, 1);
% plot(MOD_DATA);
% title('MOD DATA');
% grid on;
% subplot(2, 1, 2);
% plot(0:fs/length(MOD_DATA_FFT):fs - fs/length(MOD_DATA_FFT) ,abs(MOD_DATA_FFT));
% title('MOD DATA FFT');
% grid on;


scatterplot(DATA);

figure;
plot(MOD_DATA);
title('MOD DATA');
grid on;
% �������������

q = 1; p = 1;
for i=1:size
    MOD_DATA_Interp(q:p*L)=MOD_DATA(i);
    q = p*L + 1;
    p = p + 1;
end

MOD_DATA_Interp_FFT = fft(MOD_DATA_Interp);

figure;
subplot(2, 1, 1);
plot(MOD_DATA_Interp);
title('MOD DATA Interp');
grid on;
subplot(2, 1, 2);
plot(0:fs/length(MOD_DATA_Interp_FFT):fs - fs/length(MOD_DATA_Interp_FFT) ,abs(MOD_DATA_Interp_FFT));
title('MOD DATA Interp FFT');
grid on;

% ������� ���������� ���������
t = (0: length(MOD_DATA_Interp) - 1)/fs;
% ��������� �������� �������������� ������� � �������� fc (carrying)
car_sig = sin(2*pi*fc*t);



% ������� ������������������ MOD_DATA_interp �� ������� ������� ����� ������������ �� ������������� ������
for i = 1:length(MOD_DATA_Interp)
  SENT_TO_WAVEFORM_GENERATOR(i) = MOD_DATA_Interp(i)*car_sig(i);
end


% ��������� ������ � ������ ������� SENT_TO_WAVEFORM_GENERATOR

SENT_TO_WAVEFORM_GENERATOR_FFT = fft(SENT_TO_WAVEFORM_GENERATOR);
figure;
subplot(2, 1, 1);
plot(SENT_TO_WAVEFORM_GENERATOR);
title('SENT TO WAVEFORM GENERATOR');
grid on;
subplot(2, 1, 2);
plot(0:fs/length(SENT_TO_WAVEFORM_GENERATOR_FFT):fs - fs/length(SENT_TO_WAVEFORM_GENERATOR_FFT) ,abs(SENT_TO_WAVEFORM_GENERATOR_FFT));
title('SENT TO WAVEFORM GENERATOR FFT');
grid on;

%% ����������� � ���������� � �������� ������� �� ����

% Find a VISA-USB object.
% ������ �� �������� ����� ������ � ���� �� ������� ����������
% ���� ������ ��� ���������� �� �� ������� � ���������� obj1, �� ����������
% ��� ���� (04.10.2021, 1:42)
WG_obj = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x0957::0x2807::MY57401328::0::INSTR', 'Tag', '');

% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(WG_obj)
    WG_obj = visa('KEYSIGHT', 'USB0::0x0957::0x2807::MY57401328::0::INSTR');
else
    fclose(WG_obj);
    WG_obj = WG_obj(1);
end

name = 'my_waveforms';
sRate = fs;
amp = 1;

% Connect to instrument object, obj1.
% fopen(obj1);

% vAddress = ['USB0::0x0957::0x2807::MY57401329::0::INSTR']; %build visa address string to connect
% fgen = visa('AGILENT',vAddress); %build IO object
% obj1.Timeout = 15; %set IO time out
%calculate output buffer size
obj1_buffer = length(SENT_TO_WAVEFORM_GENERATOR)*8;
set (WG_obj,'OutputBufferSize',(obj1_buffer+125));

WG_obj.Timeout = 10;

%open connection to 33500A/B waveform generator
try
   fopen(WG_obj);
catch exception %problem occurred throw error message
    uiwait(msgbox('Error occurred trying to connect to the 33522, verify correct IP address','Error Message','error'));
    rethrow(exception);
end

%Query Idendity string and report
fprintf (WG_obj, '*IDN?');
idn = fscanf (WG_obj);
fprintf (idn)
fprintf ('\n\n')

%create waitbar for sending waveform to 33500
mes = ['Connected to ' idn ' sending waveforms.....'];
h = waitbar(0,mes);

%Reset instrument
fprintf (WG_obj, '*RST');

%make sure waveform data is in column vector
if isrow(SENT_TO_WAVEFORM_GENERATOR) == 0
    SENT_TO_WAVEFORM_GENERATOR = SENT_TO_WAVEFORM_GENERATOR';
end

%set the waveform data to single precision
SENT_TO_WAVEFORM_GENERATOR = single(SENT_TO_WAVEFORM_GENERATOR);

ON_OFF_FILTER_CH1 = ['SOURce1:FUNCtion:ARBitrary:FILTer ', 'OFF'];
fprintf(WG_obj, ON_OFF_FILTER_CH1); % ON OFF filter

%scale data between 1 and -1
mx = max(abs(SENT_TO_WAVEFORM_GENERATOR));
SENT_TO_WAVEFORM_GENERATOR = (1*SENT_TO_WAVEFORM_GENERATOR)/mx;

%update waitbar
waitbar(.1,h,mes);

%send waveform to 33500
fprintf(WG_obj, 'SOURce1:DATA:VOLatile:CLEar'); %Clear volatile memory
fprintf(WG_obj, 'FORM:BORD SWAP');  %configure the box to correctly accept the binary arb points
SENT_TO_WG_Bytes=num2str(length(SENT_TO_WAVEFORM_GENERATOR) * 4); %# of bytes
header= ['SOURce1:DATA:ARBitrary ' name ', #' num2str(length(SENT_TO_WG_Bytes)) SENT_TO_WG_Bytes]; %create header
binblockBytes = typecast(SENT_TO_WAVEFORM_GENERATOR, 'uint8');  %convert datapoints to binary before sending
fwrite(WG_obj, [header binblockBytes], 'uint8'); %combine header and datapoints then send to instrument
fprintf(WG_obj, '*WAI');   %Make sure no other commands are exectued until arb is done downloadin
%update waitbar
waitbar(.8,h,mes);
%Set desired configuration for channel 1
command = ['SOURce1:FUNCtion:ARBitrary ' name];
%fprintf(fgen,'SOURce1:FUNCtion:ARBitrary GPETE'); % set current arb waveform to defined arb testrise
fprintf(WG_obj,command); % set current arb waveform to defined arb testrise
command = ['MMEM:STOR:DATA1 "INT:\' name '.arb"'];
%fprintf(fgen,'MMEM:STOR:DATA1 "INT:\GPETE.arb"');%store arb in intermal NV memory
fprintf(WG_obj,command);
%update waitbar
waitbar(.9,h,mes);
command = ['SOURCE1:FUNCtion:ARB:SRATe ' num2str(sRate)]; %create sample rate command
fprintf(WG_obj,command);%set sample rate
fprintf(WG_obj,'SOURce1:FUNCtion ARB'); % turn on arb function
command = ['SOURCE1:VOLT ' num2str(amp)]; %create amplitude command
fprintf(WG_obj,command); %send amplitude command
fprintf(WG_obj,'SOURCE1:VOLT:OFFSET 0'); % set offset to 0 V
fprintf(WG_obj,'OUTPUT1 ON'); %Enable Output for channel 1
fprintf('SENT_TO_WG waveform downloaded to channel 1\n\n') %print waveform has been downloaded

%get rid of message box
waitbar(1,h,mes);
delete(h);

%Read Error
fprintf(WG_obj, 'SYST:ERR?');
errorstr = fscanf (WG_obj);

% error checking
if strncmp (errorstr, '+0,"No error"',13)
   errorcheck = 'Arbitrary waveform generated without any error\n';
   fprintf (errorcheck)
else
   errorcheck = ['Error reported: ', errorstr];
   fprintf (errorcheck)
end

fclose(WG_obj);

%% ����������� � ������������ � ��������� ������ �� ����

close all;
clc;
t = menu('���� ������� ������� ��������� (���)', '20', '50', '100', '200', '500', '1ms', '2ms', '5ms');


OSCI_Obj = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x2A8D::0x1797::CN58056332::0::INSTR', 'Tag', '');
% Create the VISA-USB object if it does not exist
% otherwise use the object that was found.
if isempty(OSCI_Obj)
    OSCI_Obj = visa('KEYSIGHT', 'USB0::0x2A8D::0x1797::CN58056332::0::INSTR');
else 
    fclose(OSCI_Obj);
    OSCI_Obj = OSCI_Obj(1);
end

% Set the buffer size
OSCI_Obj.InputBufferSize = 1000000;
% Set the timeout value
OSCI_Obj.Timeout = 10;
% Set the Byte order
OSCI_Obj.ByteOrder = 'littleEndian';
% Open the connection
fopen(OSCI_Obj);
% Instrument control and data retreival

% Now control the instrument using SCPI commands. refer to the instrument
% programming manual for your instrument for the correct SCPI commands for
% your instrument.

% Reset the instrument and autoscale and stop
% fprintf(visaObj,'*RST; :AUTOSCALE'); 
 fprintf(OSCI_Obj,':STOP');
% Specify data from Channel 1
fprintf(OSCI_Obj,':WAVEFORM:SOURCE CHAN1'); 
% Set timebase to main
fprintf(OSCI_Obj,':TIMEBASE:MODE MAIN');
% Set up acquisition type and count. 
fprintf(OSCI_Obj,':ACQUIRE:TYPE NORMAL');
fprintf(OSCI_Obj,':ACQUIRE:COUNT 1');
% Specify 5000 points at a time by :WAV:DATA?
fprintf(OSCI_Obj,':WAV:POINTS:MODE RAW');
fprintf(OSCI_Obj,':WAV:POINTS 50000');
% Now tell the instrument to digitize channel1
fprintf(OSCI_Obj,':DIGITIZE CHAN1');
% Wait till complete
operationComplete = str2double(query(OSCI_Obj,'*OPC?'));
while ~operationComplete
    operationComplete = str2double(query(OSCI_Obj,'*OPC?'));
end
% Get the data back as a WORD (i.e., INT16), other options are ASCII and BYTE
fprintf(OSCI_Obj,':WAVEFORM:FORMAT WORD');
% Set the byte order on the instrument as well
fprintf(OSCI_Obj,':WAVEFORM:BYTEORDER LSBFirst');
% Get the preamble block
preambleBlock = query(OSCI_Obj,':WAVEFORM:PREAMBLE?');
% The preamble block contains all of the current WAVEFORM settings.  
% It is returned in the form <preamble_block><NL> where <preamble_block> is:
%    FORMAT        : int16 - 0 = BYTE, 1 = WORD, 2 = ASCII.
%    TYPE          : int16 - 0 = NORMAL, 1 = PEAK DETECT, 2 = AVERAGE
%    POINTS        : int32 - number of data points transferred.
%    COUNT         : int32 - 1 and is always 1.
%    XINCREMENT    : float64 - time difference between data points.
%    XORIGIN       : float64 - always the first data point in memory.
%    XREFERENCE    : int32 - specifies the data point associated with
%                            x-origin.
%    YINCREMENT    : float32 - voltage diff between data points.
%    YORIGIN       : float32 - value is the voltage at center screen.
%    YREFERENCE    : int32 - specifies the data point where y-origin
%                            occurs.
% Now send commmand to read data
fprintf(OSCI_Obj,':WAV:DATA?');
% read back the BINBLOCK with the data in specified format and store it in
% the waveform structure. FREAD removes the extra terminator in the buffer
waveform.RawData = binblockread(OSCI_Obj,'uint16'); fread(OSCI_Obj,1);
% Read back the error queue on the instrument
instrumentError = query(OSCI_Obj,':SYSTEM:ERR?');
while ~isequal(instrumentError,['+0,"No error"' newline])
    disp(['Instrument Error: ' instrumentError]);
    instrumentError = query(OSCI_Obj,':SYSTEM:ERR?');
end

RECIEVED_FROM_OSCI = waveform.RawData;

%% ��������� ���������� ������ � ������������


% ������ � ������ ����������� � ������������ �������

% RECIEVED_FROM_OSCI_FFT = fft(RECIEVED_FROM_OSCI);
% figure;
% subplot(2, 1, 1);
% plot(RECIEVED_FROM_OSCI);
% title('RECIEVED FROM OSCI');
% grid on;
% subplot(2, 1, 2);
% plot(0:fs/length(RECIEVED_FROM_OSCI_FFT):fs - fs/length(RECIEVED_FROM_OSCI_FFT) ,abs(RECIEVED_FROM_OSCI_FFT));
% title('RECIEVED FROM OSCI FFT');
% grid on;



switch t
    case 1
        deci_factor_for_received_data = 5;
        deci_factor_for_sent_data = 1;
    case 2
        deci_factor_for_received_data = 2;
        deci_factor_for_sent_data = 1;
    case 3
        deci_factor_for_received_data = 1;
        deci_factor_for_sent_data = 1;
    case 4
        deci_factor_for_received_data = 1;
        deci_factor_for_sent_data = 2;
    case 5
        deci_factor_for_received_data = 1;
        deci_factor_for_sent_data = 5;
    case 6
        deci_factor_for_received_data = 1;
        deci_factor_for_sent_data = 10;
    case 7
        deci_factor_for_received_data = 1;
        deci_factor_for_sent_data = 20;
    case 8
        deci_factor_for_received_data = 1;
        deci_factor_for_sent_data = 50;
    otherwise
        % ���� �� ���� �� ��������� �� ��� ������ (�.�. ������ �������), ��
        % ���������� ���������� ���������� ���������
end



decimated_RECIEVED_FROM_OSCI = RECIEVED_FROM_OSCI(deci_factor_for_received_data:deci_factor_for_received_data:end);
decimated_RECIEVED_FROM_OSCI = decimated_RECIEVED_FROM_OSCI/max(abs(decimated_RECIEVED_FROM_OSCI))*2-1;

% ������������ ����������� (???)
% RECIEVED_FROM_OSCI = RECIEVED_FROM_OSCI(1:10000);
% local_L = 5;
% 
% local_q = 1; local_p = 1;
% for i=1:length(RECIEVED_FROM_OSCI)
%     decimated_RECIEVED_FROM_OSCI(local_q:local_p*local_L)=RECIEVED_FROM_OSCI(i);
%     local_q = local_p*local_L + 1;
%     local_p = local_p + 1;
% end


% figure(221);
% subplot(2, 1, 1);
% plot(SENT_TO_WAVEFORM_GENERATOR(700:1000));
% grid on;
% title('SENT TO WAVEFORM GENERATOR');
% subplot(2, 1, 2);
% plot(decimated_RECIEVED_FROM_OSCI);
% grid on;
% title(sprintf('decimated RECIEVED FROM OSCI with factor %f', decimation_factor));


% decimated_RECIEVED_FROM_OSCI_FFT = fft(decimated_RECIEVED_FROM_OSCI);

% figure;
% subplot(2, 1, 1);
% plot(decimated_RECIEVED_FROM_OSCI);
% title('decimated RECIEVED FROM OSCI');
% grid on;
% subplot(2, 1, 2);
% plot(0:fs/length(decimated_RECIEVED_FROM_OSCI_FFT):fs - fs/length(decimated_RECIEVED_FROM_OSCI_FFT) ,abs(decimated_RECIEVED_FROM_OSCI_FFT));
% title('decimated RECIEVED FROM OSCI FFT');
% grid on;

% ���������� ����������
SENT_TO_WAVEFORM_GENERATOR_DECIMED = SENT_TO_WAVEFORM_GENERATOR(deci_factor_for_sent_data:deci_factor_for_sent_data:end);

[c, lags] = xcorr(decimated_RECIEVED_FROM_OSCI, SENT_TO_WAVEFORM_GENERATOR_DECIMED);

rt= xcorr(decimated_RECIEVED_FROM_OSCI, SENT_TO_WAVEFORM_GENERATOR_DECIMED);
figure;
plot(abs(rt));

[~, O] = max(abs(rt));
t_lag2 = lags(O);


% ������������
c = c/max(c);
[~, I] = max(c);
t_lag = lags(I);

figure;
subplot(2, 1, 1);
plot(c);
title('CORRELATION C');
grid on;
subplot(2, 1, 2);
plot(lags);
title('CORRELATION lags');
grid on;
% 
figure;
plot(lags, abs(c));
title('CORRELATION');
grid on;

data_piece = decimated_RECIEVED_FROM_OSCI(abs(t_lag2):abs(t_lag2) + length(SENT_TO_WAVEFORM_GENERATOR_DECIMED)-1);

% figure;
% subplot(2, 1, 1);
% plot(SENT_TO_WAVEFORM_GENERATOR);
% title('SENT TO WAVEFORM GENERATOR');
% grid on;
% subplot(2, 1, 2);
% plot(data_piece);
% title('DATA PIECE');
% grid on;

% ���� ��������� ������������

% L2 = deci_factor_for_sent_data; 

% ����� ��������� ������� if deci_factor_for_sent_data > 1, ����� ��
% ��������� ����� ���� �������, ��� deci_factor_for_sent_data = 1, �� �
% ����� ������ ���������� data_piece_Interp ��� ����� �������� ������ �� ����������
% ����� ������ ������� data_piece_Interp = data_piece;

q = 1; p = 1;
for i=1:length(data_piece)
    data_piece_Interp(q:p*deci_factor_for_sent_data)=data_piece(i);
    q = p*deci_factor_for_sent_data + 1;
    p = p + 1;
end

% ������� �� ������� �������

for i = 1:length(data_piece_Interp)
  REVIVE(i) = data_piece_Interp(i)*car_sig(i);  % ������� ������� �����. �� ������� � ����� ������� ������ 5000 � ������������� L = 50,
  % � ������ �������������� ������� ������ 2500 � ������������� L = 25; ������� ����� ����
  % ������ ����� ����� ��������������� data_piece � ������������� L2 = 2
end



% ������ ������������� �� 0 ������� ��������� ������� � ��� ������
% REVIVE_FFT = fft(REVIVE);

% figure;
% subplot(2, 1, 1);
% plot(REVIVE);
% title('REVIVE');
% grid on;
% subplot(2, 1, 2);
% plot(0:fs/length(REVIVE_FFT):fs - fs/length(REVIVE_FFT) ,abs(REVIVE_FFT));
% title('RREVIVE FFT');
% grid on;

% ���������� �� ������� ������� ���� ��������������

for i = 1:length(REVIVE)/L
  a(i) = sum(REVIVE(i*L - L+1:i*L))/(L/2); 
  y(i*L - L + 1:L*i) = a(i);               
end




% ������ ���������������� ������� 
% y_FFT = fft(y);

% figure;
% subplot(2, 1, 1);
% plot(y);
% title('y');
% grid on;
% subplot(2, 1, 2);
% plot(0:fs/length(y_FFT):fs - fs/length(y_FFT) ,abs(y_FFT));
% title('y FFT');
% grid on;

% ��������� � �����������

decim = y(L:L:end);

% decim_FFT = fft(decim);

% figure;
% subplot(2, 1, 1);
% plot(decim);
% title('decim');
% grid on;
% subplot(2, 1, 2);
% plot(0:fs/length(decim_FFT):fs - fs/length(decim_FFT) ,abs(decim_FFT));
% title('decim FFT');
% grid on;

% ����������� �������������� ������ � ������� �����

for i = 1:length(decim)
    if decim(i) < 0
        demod(i) = 0;
    else
        demod(i) = 1;
    end
end

% figure(8)
% plot(decim);
% grid on
% title('������ ��������� ���������������');


% figure
% plot(demod);
% grid on;
% title('demod');
% 
% figure(111);
% subplot(2, 1, 1);
% plot(RECIEVED_FROM_OSCI);
% title('RECIEVED FROM OSCI');
% grid on;
% subplot(2, 1, 2);
% plot(SENT_TO_WAVEFORM_GENERATOR);
% title('SENT TO WAVEFORM GENERATOR');
% grid on;


scatterplot(decim);

err = biterr(DATA, demod);

[number, ratio] = biterr(DATA, demod);