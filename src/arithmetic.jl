

"""
Create tensors for MPO representing addition
The first tensor corresponds to the least significant digit.

Each tensor has two physical indices of size 2, a_i, b_i (=0, 1).

A = M^(1)_{1,d_1}(a_1, b_1) M^(2)_{d_1,d_2}(a_2, b_2) ...  M^(Q)_{d_{Q-1},d_Q}(a_Q, b_Q)
"""


"""
MPO tensor representing adder
Return size: (left bod dim, right bond dim, input a, input b, output)
"""
function _adder(;half::Bool=false, last::Bool=false, antiperiodic::Bool=true)
    cin_size = half ? 1 : 2
    cout_size = last ? 1 : 2
    adder = Tensor(ComplexF64, cin_size, cout_size, 2, 2, 2)
    for cin in 0:(cin_size-1), b in 0:1, a in 0:1
        res = a + b + cin
        cout = (res & 0b10) >> 1
        if last
            # Apply antiperiodic bc
            adder[cin+1, 1, a+1, b+1, (res & 1)+1] = (cout != 0 && antiperiodic) ? -1 : 1
        else
            adder[cin+1, cout+1, a+1, b+1, (res & 1)+1] = 1
        end
    end
    return adder
end

"""
MPO tensor selecting a or b
Return size: (left bod dim=1, right bond dim=1, input a, input b, output)
"""
function _selector(which=:left)
    (which == :left || which == :right) || error("which must be :left or :right")
    cin_size = 1
    cout_size = 1
    selector = Tensor(ComplexF64, cin_size, cout_size, 2, 2, 2)
    for b in 0:1, a in 0:1
        res = (which == :left ? a : b)
        selector[1, 1, a+1, b+1, res+1] = 1
    end
    return selector
end

function fouriertransformer()
    #
end