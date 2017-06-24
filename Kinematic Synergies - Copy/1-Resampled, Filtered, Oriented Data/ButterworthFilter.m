function [filtered] = ButterworthFilter(SF, fc, ResampledData)

% get butter worth filter coefficients
[b,a]=butter(2,fc/(SF/2));
filtered=filtfilt(b,a,ResampledData.data);

end
