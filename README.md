# int_kidney_exchange
Computing international kidney exchange schemes

This repository provides access to codes and testgames we used during our simulations in various published [[1]](#1) and still ongoing research [[2]](#2).

Below we provide the pseudo-code of the greedy search version of ${\tt Lex-Min}$ that we implemented for our computational experiments. The difference of this pseudo-code and the pseudo-code of ${\tt Lex-Min}$ in [[1]](#1) and [[3]](#3) is that in the pseudo-code outlined here we combine Steps 1 and 2 of ${\tt Lex-Min}$ together and we do the same for Steps 3 and 4 of ${\tt Lex-Min}$.

> **Input:** A partitioned matching game $(N,v)$ and an allocation $x$.<br>
> **Output:** A matching $M\in {\cal M}$ that is lexicographically minimal for $x$.<br>
> let $r=1$ \% $r-1$ *is the number of countries already finished*\;<br>
> compute an arbitrary maximum matching $M\in {\cal M}$\;<br>
> for every country $p$, let $\delta_p=|x_p-s_p(M)|$\;<br>
> let $\pi=(\pi_1,\ldots,\pi_n)$ be an order of the countries by weakly decreasing $\delta$-values\;<br>
> **while** $r \leq n$ \%*for the unfinished countries* **do**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;**if** $\delta_{\pi_r}\leq 0.5$ **then**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;STOP and RETURN $M$ \% $M$ *is lexicographically minimal*\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;**endif**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;let $R_{\pi_r}=\ ]x_{\pi_r}-\delta_{\pi_r}, x_r+\delta_{\pi_r}[$ \% *trying to reduce the difference of country* ${\pi_r}$\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;**for** $k=r+1,\dots, n$ **do**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**if** $\delta_{\pi_k}=\delta_{\pi_r}$ **then**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;let $R_{\pi_k}=[x_{\pi_k}-\delta_{\pi_r}, x_{\pi_k}+\delta_{\pi_r}]$\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**else**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;let $R_{\pi_k}=\ ]x_{\pi_k}-\delta_{\pi_r}, x_{\pi_k}+\delta_{\pi_r}[$\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**endif**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;**endfor**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;let $I_{\pi_k}$ be the largest non-negative integer sub-interval of $R_{\pi_k}, k=r, \dots, n$\; <br>
> &nbsp;&nbsp;&nbsp;&nbsp;**if** Corollary 3.1 of [[1]](#1) returns a matching $M'$ for intervals $I_{\pi_r},\ldots,I_{\pi_n}$ **then**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;let $M=M'$ and recompute $\pi$\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\% $\pi$ *is a new order with respect to the new deviations* $|x_p-s_p(M')|$\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;**else**<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\% $\delta_{\pi_r}$ *is fixed and country* ${\pi_r}$ *is finished, we keep* $\pi$ *and* $M$ *the same*\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;let $r=r+1$\;<br>
> &nbsp;&nbsp;&nbsp;&nbsp;**endif**<br>
> **endwhile**<br>

## References
<a id="1">[1]</a> 
M. Benedek, P. Biró, W. Kern and D. Paulusma (2022).
[Computing Balanced Solutions for Large International Kidney Exchange Schemes.](https://dl.acm.org/doi/10.5555/3535850.3535861)<br>
*AAMAS '22: Proceedings of the 21st International Conference on Autonomous Agents and Multiagent Systems*, 82&ndash;90.

<a id="2">[2]</a> 
M. Benedek, P. Biró, D. Paulusma and X. Ye (2023).
[Computing Balanced Solutions for Large International Kidney Exchange Schemes.](https://arxiv.org/abs/2109.06788)<br>
*arXiv preprint 2109.06788*.

<a id="3">[3]</a> 
M. Benedek, P. Biró, W. Kern, D. Pálvölgyi and D. Paulusma (2023).
[Partitioned Matching Games for International Kidney Exchange.](https://arxiv.org/abs/2301.13181)<br>
*arXiv preprint 2301.13181*.
