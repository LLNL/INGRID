import geometry as g
import matplotlib.pyplot as plt



p1 = g.Point((1.235e+00, 1.215e+00))
p2 = g.Point((1.235e+00, 1.205e+00))
L1 = g.Line([p1, p2])

q1 = g.Point((1.234e+00, 1.211e+00))
q2 = g.Point((1.236e+00, 1.21e+00))
L2 = g.Line([q1, q2])

plt.figure()
L1.plot('blue')
L2.plot('red')

res, sol = g.segment_intersect(L1.points(), L2.points())
plt.plot(sol[0][0], sol[0][1], 'x', color = 'cyan', zorder = 10)
plt.plot(sol[1][0], sol[1][1], '*', color = 'green', zorder = 10)

plt.show()