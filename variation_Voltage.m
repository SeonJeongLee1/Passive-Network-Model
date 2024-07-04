function [Vset_out, Vreset_out] = variation_Voltage(Vset, Vreset, Vset_change, Vreset_change, num)
    Vset_out = zeros(num, 1);
    Vreset_out = zeros(num, 1);
    for i = 1:num
        Vset_out(i) = Vset + Vset_change*randn(1,1);
        Vreset_out(i) = Vreset + Vreset_change*randn(1,1);
    end
end