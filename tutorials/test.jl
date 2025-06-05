using Plots

p = plot(1, xlim=(0, 3π), ylim=(-1.5, 1.5), title="Sine", marker=2)
N = 100
@gif for i=1:N
    x = (i-1) * 3π / N
    push!(p, x, sin(x))
end