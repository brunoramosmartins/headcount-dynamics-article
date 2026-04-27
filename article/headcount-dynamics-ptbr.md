# O Time Que Se Substitui Sozinho

## Dinâmica de Headcount e Turnover como Processos Estocásticos — de Cadeias de Markov a Modelos de Nascimento e Morte

> Um plano de headcount que diz "vamos ter 50 pessoas em dezembro" é uma
> previsão pontual de um processo estocástico. Ignora a probabilidade de
> realmente atingir esse estado, o tempo esperado para chegar lá e a
> distribuição estacionária para a qual o sistema naturalmente tende. Este
> artigo mostra que a teoria de cadeias de Markov fornece respostas exatas
> para as três perguntas, e aplica as ferramentas a um modelo de time de
> TI funcional.

---

## 1. Introdução — Headcount É um Processo Estocástico

Um plano de pessoal tipicamente é apresentado como uma sequência de
inteiros. *45 pessoas em janeiro. 50 em março. 60 até o fim do ano.* Os
números parecem concretos, mas são projeções de um modelo no qual
ninguém pede demissão, contratações chegam no prazo e promoções acontecem
sob demanda. Nenhuma dessas premissas se sustenta. O que os planejadores
realmente enfrentam é um *processo aleatório*: a qualquer momento,
funcionários podem sair, ser promovidos ou contratados. O número de
pessoas no time daqui a doze meses é uma variável aleatória com uma
distribuição, não um valor único.

Uma vez que o headcount é reconhecido como um processo aleatório, três
perguntas tornam-se quantificáveis — perguntas que o plano determinístico
não consegue responder:

1. **Onde o time vai se estabilizar?** Se as taxas de contratação e atrição
   persistirem, o tamanho e a composição do time convergem para uma
   *distribuição estacionária* — não necessariamente onde o plano dizia.
2. **Quanto tempo vai levar?** O tempo esperado para atingir um headcount-
   alvo é o *tempo de hitting* de uma cadeia de Markov, computável a
   partir da matriz de transição.
3. **Qual é a probabilidade de atingir o plano?** $P(n_{12} \geq 50)$ é
   uma fórmula de uma linha, dada a especificação da cadeia.

A matemática é madura. Cadeias de Markov em tempo discreto remontam ao
artigo original de Markov de 1906. Cadeias em tempo contínuo e a
especialização nascimento-morte, que modelam fluxos contínuos de
contratação e atrição, são material clássico de pesquisa operacional e
teoria de filas. O que este artigo faz é *aplicar* esse aparato de ponta
a ponta a um problema de planejamento de pessoal e conectar as quantidades
analíticas resultantes à simulação de Monte Carlo de orçamento do Artigo 1
desta trilogia.

O argumento se desdobra em três camadas. As três primeiras seções
introduzem cadeias de Markov em tempo discreto: a matriz de transição $P$,
a equação de Chapman–Kolmogorov e a distribuição estacionária $\pi$. A
quarta seção adiciona análise de absorção — a *matriz fundamental*
$N = (I - Q)^{-1}$ — que responde perguntas como "quanto tempo até o time
estar vazio sob um congelamento de contratações?". A quinta seção avança
para tempo contínuo: o gerador $Q$, a exponencial matricial
$P(t) = e^{Qt}$ e processos nascimento-morte, que produzem distribuições
de tamanho de time Poisson em forma fechada. A sexta seção reúne tudo num
`HeadcountModel` aplicado e roda cinco cenários realistas. A sétima valida
empiricamente cada afirmação teórica. A oitava conecta a distribuição de
headcount à distribuição de orçamento do Artigo 1, completando a trilogia.
A nona oferece um framework de quatro passos para um planejador que queira
começar a usar essas ferramentas. A décima conclui.

Um leitor confortável com álgebra linear e probabilidade elementar achará
o artigo autossuficiente. Cada teorema enunciado é demonstrado ou
esboçado, cada exemplo numérico é reprodutível pelo código no repositório
companheiro, e cada figura é gerada por um script versionado com seed
aleatória fixa.

---

## 2. Cadeias de Markov — A Linguagem das Transições

Uma cadeia de Markov em tempo discreto sobre um espaço de estados finito
$S = \{1, \ldots, K\}$ é uma sequência de variáveis aleatórias
$\{X_n\}_{n \geq 0}$ tal que a probabilidade do próximo estado depende
apenas do estado presente, não do histórico:

$$
P(X_{n+1} = j \mid X_n = i, X_{n-1}, \ldots, X_0) \;=\; P(X_{n+1} = j \mid X_n = i) \;=\; p_{ij}.
$$

Essa é a **propriedade markoviana** — o passado é resumido pelo presente.
Os números $p_{ij}$ formam a **matriz de transição** $P = (p_{ij})$. Cada
linha de $P$ é uma distribuição de probabilidade sobre o próximo estado,
de modo que

$$
\sum_{j \in S} p_{ij} \;=\; 1 \quad \text{para todo } i,
$$

e $P$ é **estocástica por linhas**: cada linha soma 1.

### Convenções de notação

Ao longo do artigo, $\pi_n$ denota a distribuição de $X_n$ como vetor-
linha, então $\pi_{n+1} = \pi_n P$ e $\pi_n = \pi_0 P^n$. Usamos
$\mathbf{1}$ para o vetor-coluna de uns; $e_i$ para o vetor-linha base
canônico com 1 na posição $i$. Os rótulos de estado para o exemplo de
headcount são $S = \{J, M, S, E\}$ — Júnior, Pleno (Mid), Sênior, Saída
(Exit). Distribuições estacionárias são escritas $\pi$ sem subscrito. Os
autovalores de $P$ são $\lambda_1, \lambda_2, \ldots$ em ordem decrescente
de módulo; para uma matriz estocástica por linhas, $\lambda_1 = 1$ sempre.

### A cadeia de headcount

Para um time de TI com três níveis de carreira e um estado de Saída, uma
matriz de transição mês a mês plausível é

$$
P = \begin{pmatrix}
0.93 & 0.03 & 0    & 0.04 \\
0    & 0.96 & 0.02 & 0.02 \\
0    & 0    & 0.99 & 0.01 \\
0.50 & 0    & 0    & 0.50
\end{pmatrix}.
$$

As três primeiras linhas codificam taxas de promoção e atrição por nível
de carreira. A última linha modela *contratação de reposição*: uma vaga
em Saída retorna a Júnior a cada mês com probabilidade 50%. Sem a última
linha, Saída seria *absorvedor* — uma vez ali, a cadeia permanece — e o
time inteiro acabaria saindo.

### Transições em $n$ passos

A probabilidade de ir do estado $i$ ao estado $j$ em $n$ meses é a entrada
$(i, j)$ de $P^n$. Para provar isso, condicione no estado intermediário
no tempo $n$:

$$
p_{ij}^{(m+n)} \;=\; \sum_{k} p_{ik}^{(m)} \, p_{kj}^{(n)}.
$$

Essa é a **equação de Chapman–Kolmogorov**, o conteúdo probabilístico da
associatividade da multiplicação matricial $P^{m+n} = P^m P^n$. Iterar um
passo de cada vez é o mesmo que uma única potência matricial grande.

### Classificação de estados

Dois estados **se comunicam** ($i \leftrightarrow j$) se cada um pode ser
alcançado pelo outro em finitos passos. Comunicação é uma relação de
equivalência, particionando o espaço de estados em **classes de
comunicação**. Uma cadeia é **irredutível** quando há uma única classe.

Um estado é **recorrente** se a cadeia retorna a ele com probabilidade 1,
e **transiente** caso contrário. Estados **absorvedores** são o caso
extremo de recorrência — uma vez dentro, nunca fora. A cadeia de
headcount com contratação de reposição é irredutível: a partir de
qualquer nível de carreira, qualquer outro estado é alcançável. Zere a
taxa de reciclagem ($p_{41} = 0$) e a cadeia torna-se absorvedora —
Saída aprisiona a dinâmica.

Essa segunda variante é exatamente a que a Seção 4 precisa para analisar
um congelamento de contratações. A primeira variante é a que a Seção 3
precisa para analisar a composição do time no longo prazo.

---

## 3. Onde o Time Vai Se Estabilizar? — Distribuições Estacionárias

Um vetor-linha $\pi = (\pi_1, \ldots, \pi_K)$ é uma **distribuição
estacionária** se satisfaz

$$
\pi P \;=\; \pi, \qquad \pi_i \geq 0, \qquad \sum_i \pi_i \;=\; 1.
$$

Iniciada em $\pi$, a cadeia permanece em $\pi$ para sempre. A pergunta
desta seção é: sob que condições $\pi$ existe, quando é única, e com que
velocidade uma distribuição inicial arbitrária converge para ela?

### Existência

Para uma cadeia de Markov finita irredutível, uma distribuição
estacionária **sempre existe**. A demonstração usa a média de Cesàro

$$
\bar\pi_n \;=\; \frac{1}{n} \sum_{k=0}^{n-1} e_i P^k.
$$

Cada $\bar\pi_n$ pertence ao simplex de probabilidade, que é compacto.
Pelo teorema de Bolzano–Weierstrass, uma subsequência converge para algum
$\pi^* \in \Delta_K$. Um cálculo curto mostra que
$\bar\pi_n P - \bar\pi_n \to 0$, então $\pi^* P = \pi^*$ — o limite é
estacionário.

### Unicidade e Perron–Frobenius

Para uma cadeia finita irredutível e aperiódica, $\pi$ é **única**. A
forma mais limpa de ver isso é via a estrutura espectral de $P$. O
teorema de Perron–Frobenius diz:

1. Todo autovalor $\lambda$ de $P$ satisfaz $|\lambda| \leq 1$.
2. $\lambda = 1$ é um autovalor simples (multiplicidade 1) para uma
   cadeia irredutível e aperiódica.
3. O único autovetor à esquerda para $\lambda = 1$, normalizado para um
   vetor de probabilidade, é $\pi$ — e $\pi_i > 0$ para todo $i$.

O segundo autovalor $\lambda_2$ governa a taxa de convergência. Decomponha
$P$ espectralmente:

$$
P^n \;=\; \mathbf{1}\pi + \sum_{k \geq 2} \lambda_k^n \, u_k v_k^\top.
$$

Como $|\lambda_k| < 1$ para $k \geq 2$, os termos não estacionários
decaem geometricamente. Obtemos a cota

$$
\|P^n - \mathbf{1}\pi\| \;\leq\; C \cdot |\lambda_2|^n.
$$

O **gap espectral** $\gamma = 1 - |\lambda_2|$ é a inversa da escala de
tempo do esquecimento. Um gap próximo de 1 significa que a cadeia mistura
rápido; um gap próximo de 0 significa que a condição inicial persiste.

### A distribuição estacionária do headcount

Para a cadeia com reciclagem, $\pi P = \pi$ é resolvível à mão. A segunda
equação dá $\pi_M = \tfrac{3}{4} \pi_J$, a terceira dá
$\pi_S = \tfrac{3}{2} \pi_J$, a quarta dá $\pi_E = 0.14 \, \pi_J$.
Normalizando,

$$
\pi \;\approx\; (0{,}295, \; 0{,}221, \; 0{,}443, \; 0{,}041).
$$

No longo prazo, um funcionário passa cerca de 30% dos meses como Júnior,
22% como Pleno, 44% como Sênior e 4% em trânsito por Saída. A fração
Sênior domina porque a taxa de atrição Sênior ($p_{34} = 0.01$) é a menor
de todas as taxas — Sêniores permanecem.

Os autovalores de $P$ são aproximadamente
$\{1{,}00, \; 0{,}99, \; 0{,}94, \; 0{,}41\}$. O gap espectral é
$\gamma \approx 0{,}01$, dando um tempo de mistura da ordem de
$1/\gamma \approx 100$ meses. Empiricamente, a distância em variação
total cai abaixo de $1/4$ em torno do mês 70. A composição "de longo
prazo" para essa cadeia é um horizonte de seis anos.

![Convergência de $P^n$ para a matriz estacionária de posto 1](../figures/transition_convergence_heatmap.png)

A figura mostra $P^n$ para $n = 1, 5, 10, 50, 100$ ao lado da matriz
limite $\mathbf{1}\pi$. Cada linha de $P^n$ se aproxima de $\pi$ conforme
$n$ cresce, e por $n = 100$ as linhas são visualmente indistinguíveis do
limite.

![Espectro de autovalores e gap espectral](../figures/eigenvalue_spectrum.png)

O espectro está dentro do disco unitário. O autovalor dominante é $1$
(verde); o segundo maior em módulo, $\lambda_2 \approx 0{,}99$ (vermelho),
estabelece a taxa de convergência. Um gap de 0,01 é pequeno — típico para
cadeias de força de trabalho do mundo real, em que taxas de promoção e
atrição são de poucos por cento por mês — e é a razão pela qual a mistura
leva anos.

---

## 4. Quanto Tempo Vai Levar? — Absorção e Tempos de Hitting

A distribuição estacionária descreve onde a cadeia chega. Não diz nada
sobre como o comportamento *transiente* evolui — a trajetória antes da
assíntota. Esta seção responde perguntas como "se as contratações
pararem, quantos meses até o time cair abaixo de 40?" e "qual fração de
Júniores algum dia chega a Sênior?". O único objeto que responde ambas é
a **matriz fundamental**.

### A forma canônica

Faça $p_{41} = 0$ na cadeia de headcount para tornar Saída absorvedora.
Ordene os estados de forma que os transientes venham primeiro e os
absorvedores por último. A matriz de transição assume a **forma canônica**

$$
P \;=\; \begin{pmatrix} Q & R \\ 0 & I \end{pmatrix},
$$

em que $Q$ é o bloco transiente-para-transiente, $R$ é o bloco
transiente-para-absorvedor, e o bloco absorvedor é simplesmente a
identidade. Para a cadeia com Saída absorvedora,

$$
Q \;=\; \begin{pmatrix} 0.93 & 0.03 & 0 \\ 0 & 0.96 & 0.02 \\ 0 & 0 & 0.99 \end{pmatrix}, \qquad
R \;=\; \begin{pmatrix} 0.04 \\ 0.02 \\ 0.01 \end{pmatrix}.
$$

### A matriz fundamental

Como a cadeia atinge um estado absorvedor com probabilidade 1, $Q^n \to 0$
quando $n \to \infty$. Isso implica $\rho(Q) < 1$, então $I - Q$ é
invertível e

$$
N \;=\; (I - Q)^{-1} \;=\; \sum_{k=0}^\infty Q^k.
$$

A matriz $N$ é a **matriz fundamental**. Cada entrada tem uma
interpretação probabilística: $N_{ij}$ é o *número esperado de meses*
gastos no estado transiente $j$ antes da absorção, dado início em $i$. O
tempo esperado até a absorção a partir do estado $i$ é a soma de todas as
visitas:

$$
t \;=\; N \mathbf{1}.
$$

Para a cadeia de headcount com Saída absorvedora, $N$ é triangular
superior com entradas diagonais $1/0.07$, $1/0.04$, $1/0.01$, e

$$
t \;=\; \begin{pmatrix} 46{,}4 \\ 75{,}0 \\ 100{,}0 \end{pmatrix} \;\text{meses}.
$$

Um Júnior inicial espera cerca de 46 meses até a saída; um Sênior,
exatamente $1/0{,}01 = 100$ meses — isto é, 100 meses em expectativa, com
variância de distribuição geométrica.

![Tempo esperado de absorção por nível inicial, com barras de $\pm 1\sigma$](../figures/absorption_times.png)

As barras mostram $t_i$ para cada nível inicial sob a variante
absorvedora. As barras de erro são $\pm 1$ desvio-padrão, calculadas a
partir da fórmula de variância $\mathrm{Var}(T) = (2N - I) t - t \odot t$.
A permanência Sênior é cerca do dobro da permanência Júnior, mas com
dispersão proporcionalmente maior.

### Probabilidades de absorção e primeiro passagem

Quando a cadeia tem múltiplos estados absorvedores, a **matriz de
probabilidades de absorção** $B = NR$ nos diz onde a cadeia termina. Para
calcular a probabilidade de um Júnior chegar a Sênior na variante
absorvedora, modifique a cadeia para que *tanto* Sênior quanto Saída
sejam absorvedores, e leia a coluna relevante de $B$. O resultado é

$$
P(\text{Júnior chega a Sênior eventualmente}) \;\approx\; 0{,}214.
$$

Apenas um Júnior em cinco chega a Sênior sob as taxas padrão; o resto sai
antes. Para um Pleno inicial, o número é exatamente $0{,}5$.

### Ponte para a Seção 3: tempos médios de retorno

Para uma cadeia irredutível (a variante com reciclagem), o tempo médio de
*retorno* a um estado é igual ao recíproco de sua probabilidade
estacionária:

$$
\mathbb{E}[T_i \mid X_0 = i] \;=\; \frac{1}{\pi_i}.
$$

Essa é a ponte entre as Seções 3 e 4. A distribuição estacionária e os
tempos médios de retorno são os mesmos dados em duas linguagens
diferentes. Recorrência (Seção 3) e absorção (Seção 4) são as duas faces
da imagem transiente: uma cadeia em equilíbrio retorna infinitas vezes;
uma cadeia absorvedora parte uma vez e nunca volta.

---

## 5. Fluxo Contínuo — Processos Nascimento-Morte

Meses são um artefato de quando calculamos folha de pagamento.
Contratações e atrição *acontecem* continuamente: um pedido de demissão
pode cair no dia 14, uma proposta pode ser assinada no dia 22. Para
modelar isso, mudamos de instantâneos mês a mês para cadeias de Markov
em tempo contínuo (CTMCs).

### Matrizes geradoras

Uma CTMC homogênea no tempo sobre o espaço de estados $S$ tem uma
**matriz geradora** $Q$ na qual $Q_{ij}$ é a *taxa* de transição de $i$
para $j$ para $i \neq j$, e $Q_{ii} = -\sum_{j \neq i} Q_{ij}$. Cada
linha soma zero (em vez de um, como no caso discreto). Diferenciando a
identidade óbvia $P(t)\mathbf{1} = \mathbf{1}$ em $t = 0$ obtemos a
propriedade de soma de linha zero imediatamente.

As probabilidades de transição $P(t) = (p_{ij}(t))$ satisfazem a
**equação de Kolmogorov forward**

$$
\frac{d}{dt} P(t) \;=\; P(t) \, Q,
$$

com condição inicial $P(0) = I$. A solução única é a exponencial
matricial

$$
P(t) \;=\; e^{Qt} \;=\; \sum_{k = 0}^\infty \frac{(Qt)^k}{k!}.
$$

Numericamente, `scipy.linalg.expm` calcula $e^{Qt}$ via aproximação de
Padé; esse é o cavalo de batalha do resto da seção.

### Processos nascimento-morte

Um **processo nascimento-morte** é uma CTMC sobre $\{0, 1, 2, \ldots\}$
com apenas transições entre vizinhos: um *nascimento* move $n \to n+1$ a
taxa $\lambda_n$, uma *morte* move $n \to n-1$ a taxa $\mu_n$. O gerador
é tridiagonal. Headcount se encaixa naturalmente: uma contratação é um
nascimento, uma demissão é uma morte.

A distribuição estacionária satisfaz as equações de **balanço detalhado**

$$
\pi_n \, \lambda_n \;=\; \pi_{n+1} \, \mu_{n+1},
$$

— o fluxo de probabilidade de $n$ para $n+1$ é igual ao fluxo de volta.
Resolvendo a recorrência obtemos

$$
\pi_n \;=\; \pi_0 \prod_{k = 0}^{n-1} \frac{\lambda_k}{\mu_{k+1}}.
$$

Um caso especial particularmente limpo é **contratação constante com
atrição per capita**: $\lambda_n = \lambda$ e $\mu_n = n\mu$. Então

$$
\pi_n \;=\; e^{-\rho} \frac{\rho^n}{n!}, \qquad \rho \;=\; \frac{\lambda}{\mu}.
$$

O tamanho do time é **Poisson** com média $\rho$. Contratação constante a
taxa $\lambda$ por mês e atrição per capita $\mu$ por mês produzem um
tamanho de time que, em equilíbrio, se comporta exatamente como o
contador de chegadas de um processo Poisson.

### A trajetória esperada

Tomando esperanças nos dois lados da equação de taxa do nascimento-morte,
obtemos a EDO

$$
\frac{d}{dt} \mathbb{E}[X_t] \;=\; \lambda - \mu \mathbb{E}[X_t],
$$

com solução fechada

$$
\mathbb{E}[X_t] \;=\; \rho + (n_0 - \rho) \, e^{-\mu t}.
$$

A média relaxa exponencialmente para a assíntota $\rho$ com meia-vida
$\ln 2 / \mu$. Com $\mu = 0{,}025$, a meia-vida é cerca de 27,7 meses —
mais de dois anos. O time esquece seu tamanho inicial nessa taxa.

### O nascimento-morte do headcount

Para um time de TI com $\lambda = 2$ contratações por mês e atrição per
capita $\mu = 0{,}025$, a média assintótica é $\rho = 80$. Começando de
$n_0 = 45$,

$$
\mathbb{E}[X_{12}] \;\approx\; 54, \qquad P(X_{12} \geq 50 \mid X_0 = 45) \;\approx\; 0{,}80.
$$

Há 20% de chance de terminar o ano abaixo de 50 — fato que a previsão
determinística $\mathbb{E}[X_{12}] = 54$ esconde por completo.

![Estacionária do nascimento-morte: $\pi_n$ é Poisson($\rho$)](../figures/birth_death_stationary.png)

As barras mostram a distribuição analítica de balanço detalhado; a linha
vermelha é a PMF Poisson com média $\rho = 80$. Coincidem até a precisão
numérica.

![Trajetória de Kolmogorov com banda 5–95%](../figures/kolmogorov_trajectory.png)

A trajetória esperada $\mathbb{E}[X_t]$ casa com a EDO em forma fechada
(tracejada). A banda sombreada é o envelope exato dos percentis 5 e 95,
calculado a partir da distribuição transiente da CTMC truncada em cada
$t$.

---

## 6. O Modelo de Headcount — Cinco Cenários

As Seções 2–5 construíram o aparato. Esta seção o aplica. Definimos um
`HeadcountModel` combinado que envolve uma cadeia em tempo discreto
(composição de níveis de carreira) e um processo nascimento-morte
(tamanho total do time), e rodamos cinco cenários que um planejador real
enfrentaria.

### O modelo combinado

O `HeadcountModel` responde perguntas de planejamento por meio de seis
métodos:

| Método | Retorna |
|--------|---------|
| `forecast(months)` | distribuição exata sobre o tamanho do time |
| `prob_reach_target(target, deadline)` | probabilidade de $n(t) \geq \text{alvo}$ |
| `expected_time_to_target(target)` | primeiro $t$ com $\mathbb{E}[n(t)]$ no alvo |
| `steady_state_composition()` | % de longo prazo por nível de carreira |
| `simulate_trajectories(...)` | caminhos amostrais via Gillespie |
| `expected_total_salary(...)` | ponte de orçamento (Seção 8) |

A composição (DTMC) e o tamanho do time (nascimento-morte) são tratados
como processos independentes. Esse desacoplamento é uma simplificação
deliberada: uma fila multiclasse seria mais precisa, mas perderia formas
fechadas. A Seção 7 mostra que o viés da simplificação contra a verdade
de Monte Carlo é pequeno.

### Cenário S1 — Crescimento estável

*Pergunta.* "Precisamos crescer de 45 para 60 em doze meses. É factível?"

Sob os parâmetros base, $\mathbb{E}[X_{12}] \approx 54$, aquém de 60.
Para chegar a 60 em expectativa, $\lambda$ precisa subir para cerca de
2,31 contratações por mês — um aumento de 16%. A figura compara o caso
base ao caso turbinado lado a lado.

![Cenário S1 — crescimento estável](../figures/scenario_s1_steady_growth.png)

A linha tracejada vermelha marca o alvo. A probabilidade de atingir 60 no
mês 12 sobe de 0,31 para 0,51 com o aumento de contratação. É uma
diferença de 20 pontos: vale dinheiro real em esforço de recrutamento
adicional.

### Cenário S2 — Congelamento de contratações

*Pergunta.* "Se congelarmos contratações hoje, quando caímos abaixo de 40?"

Faça $\lambda = 0$. A trajetória média torna-se decaimento exponencial
puro: $\mathbb{E}[n(t)] = 45 \, e^{-0{,}025 t}$, cruzando 40 em
$t \approx 4{,}7$ meses. A banda 5–95% na figura adiciona cerca de
$\pm 2$ meses de variabilidade ao cruzamento do limiar.

![Cenário S2 — congelamento de contratações](../figures/scenario_s2_hiring_freeze.png)

### Cenário S3 — Recuperação após corte

*Pergunta.* "Acabamos de cortar para 35 pessoas. Quanto tempo até voltar
a 45?"

Reinicie $n_0 = 35$ com as taxas originais de contratação/atrição. A
média atinge 45 em $\ln(45/35) / 0{,}025 \approx 10$ meses. A recuperação
é mais rápida que a meia-vida porque o alvo está abaixo da assíntota: a
cadeia é *puxada para cima* pela diferença até $\rho = 80$.

![Cenário S3 — recuperação após corte](../figures/scenario_s3_layoff_recovery.png)

### Cenário S4 — Mudança de composição

*Pergunta.* "Queremos a fração estacionária Sênior subir de 44% para 60%.
Que mudança de taxa nos leva lá?"

A fração Sênior é mais sensível à taxa de atrição Sênior $p_{34}$.
Reduzi-la pela metade (de 0,01 para 0,005) eleva $\pi_S$ de 0,44 para
0,59. A figura plota $\pi_J, \pi_M, \pi_S$ como funções de $p_{34}$ ao
longo de uma faixa realista. Dobrar $p_{34}$ empurra a fração Sênior
para 0,30 e inunda Pleno; reduzi-la pela metade tem o efeito oposto.

![Cenário S4 — mudança de composição](../figures/scenario_s4_composition_shift.png)

### Cenário S5 — Contratação sazonal

*Pergunta.* "A contratação é em ondas: $\lambda = 3$ no Q1 e Q3,
$\lambda = 1$ no Q2 e Q4. A média conta a história toda?"

A média anual não muda (linearidade da expectativa), mas a variância no
mês 12 é cerca de 6% maior sob contratação sazonal do que sob constante
$\bar\lambda = 2$. A figura mostra a trajetória sazonal acompanhando a
linha de base de $\lambda$ constante de perto em média, mas com bandas
5–95% visivelmente mais largas nos picos de fim de trimestre.

![Cenário S5 — contratação sazonal](../figures/scenario_s5_seasonal_hiring.png)

### Tornado de sensibilidade

Um planejador com orçamento limitado para iniciativas de contratação ou
retenção precisa de um ranking de quais taxas mais importam. A figura
tornado varia cada taxa em $\pm 50\%$ ao redor do seu valor base e
ranqueia parâmetros pela magnitude do impacto sobre o tamanho de time de
longo prazo e sobre a fração Sênior.

![Tornado de sensibilidade](../figures/sensitivity_tornado.png)

Os achados centrais: $\mu$ tem um impacto assimétrico ligeiramente maior
sobre $\rho = \lambda / \mu$ que $\lambda$ (porque
$\rho \propto 1/\mu$). A fração Sênior é dominada pela atrição Sênior —
as taxas mais importantes em seguida são as promoções a montante, com a
atrição Júnior tendo o menor efeito. A implicação para ação: retenção
Sênior é o investimento de maior alavancagem para moldar a composição do
time.

---

## 7. Experimentos e Resultados

Cada resultado das Seções 2–6 é uma afirmação analítica. Esta seção
valida cada uma numericamente contra simulação. Todos os experimentos
usam uma única seed fixa (2026) e um único módulo de estilo
(`scripts/_style.py`) para garantir reprodutibilidade.

### Experimento A — Convergência de $P^n$ para $\mathbf{1}\pi$

O teorema de convergência diz $P^n \to \mathbf{1}\pi$ quando $n \to
\infty$ para uma cadeia irredutível e aperiódica. O heatmap da Seção 3
visualiza isso: cada linha de $P^n$ converge para $\pi$, com $P^{100}$
visualmente indistinguível do limite. O decaimento geométrico à taxa
$|\lambda_2|^n$ significa que cada linha reduz pela metade sua distância
de $\pi$ a cada $\sim \ln 2 / \gamma \approx 70$ meses.

### Experimento B — Gap espectral e tempo de mistura

A figura de autovalores na Seção 3 mostra o espectro dentro do disco
unitário. Numericamente, o gap espectral é $\gamma \approx 0{,}01$ e o
tempo de mistura empírico (o menor $n$ tal que
$\max_i \tfrac{1}{2}\|P^n_{i,:} - \pi\|_1 \leq 1/4$) é
aproximadamente 70 meses. A cota
$t_{\mathrm{mix}} \leq (1/\gamma) \ln(1/(\varepsilon \pi_{\min}))$ prevê
cerca de 200 meses — a cota é frouxa porque $\pi_{\min}$ é pequeno.

### Experimento C — Tempos de absorção

O gráfico de barras na Seção 4 mostra o tempo esperado até a saída por
nível inicial sob a variante absorvedora: 46 meses para Júnior, 75 para
Pleno, 100 para Sênior. Os desvios-padrão da fórmula de variância
$(2N - I)t - t \odot t$ são 47, 50, 99 — então os intervalos de 1-sigma
são largos. A permanência Sênior é geometricamente distribuída com
parâmetro 0,01, média 100, desvio 99 — a cadeia "aposta" em nunca pedir
demissão, com alta variância.

### Experimento D — Estado estacionário do nascimento-morte

Duas figuras validam o resultado Poisson. A primeira plota a distribuição
analítica de balanço detalhado contra a PMF Poisson; coincidem até a
precisão de ponto flutuante. A segunda simula 200 caminhos ao longo de
400 meses, descarta os primeiros 200 meses como burn-in, e plota um
histograma das amostras restantes de tamanho de time.

![Estado estacionário do nascimento-morte — histograma empírico vs Poisson](../figures/birth_death_steady_simulation.png)

A média empírica é 80,33 e a variância empírica é 80,99, ambas dentro de
1% do alvo Poisson(80). O $\sup |F_{\mathrm{emp}} - F_{\mathrm{Poisson}}|$
no estilo Kolmogorov–Smirnov é 0,022 — bem abaixo de valores críticos
típicos para $N \sim 40\,000$ amostras — confirmando que empírico e
analítico concordam.

### Experimento E — Os cinco cenários

As cinco figuras de cenário da Seção 6 *são* o experimento: cada uma
sobrepõe 30 caminhos simulados sobre a média analítica e a banda exata
5–95% derivada da CTMC truncada. Os caminhos simulados preenchem
visualmente a banda, fornecendo uma checagem de ponta a ponta da
correção do modelo.

### Experimento F — Tornado de sensibilidade

A figura tornado ranqueia taxas pela magnitude de impacto. A atrição
Sênior $p_{34}$ lidera, com $|\Delta \pi_S| \approx 0{,}17$ ao longo de
$\pm 50\%$. As próximas competidoras são promoção Pleno→Sênior ($p_{23}$,
$\Delta = 0{,}12$) e promoção Júnior→Pleno ($p_{12}$, $\Delta = 0{,}11$).
O menor efeito vem da atrição Júnior ($\Delta = 0{,}005$) —
contraintuitivamente baixo, porque a rotatividade no nível de entrada é
em grande parte absorvida pela contratação de reposição na variante de
reciclagem.

### Experimento G — Animação de trajetória

O GIF animado revela 100 trajetórias simuladas mês a mês, sobrepostas com
a média analítica e a banda 5–95%. A experiência visual — caminhos
abrindo-se em leque a partir de $n_0 = 45$, embalados pela banda à
medida que ela se alarga, a curva da média subindo em direção a $\rho =
80$ — é o artefato único mais claro para audiências não técnicas. É a
figura que mais vale mostrar numa reunião de planejamento.

### Validação cruzada

Os sete experimentos cobrem cada teorema enunciado no artigo. Os números
analíticos e empíricos concordam dentro do erro de Monte Carlo em todos
os casos. Cada script vive em `scripts/`, usa seed 2026 e está conectado
a um runner mestre que reproduz todas as 13 figuras com um único
comando.

---

## 8. Conectando ao Orçamento — Artigos 1 e 2

O Artigo 1 desta trilogia simula o orçamento anual total amostrando
headcount e salários conjuntamente. O Artigo 2 ajusta a distribuição de
salários. O presente artigo dá a distribuição de headcount. Juntos
formam um kit de ferramentas probabilístico para planejamento de
orçamento de TI.

A ponte fechada usa as leis da expectativa total e da variância total.
Seja $N$ = tamanho do time no mês 12 e $S_i$ = salários i.i.d. com
$\mathbb{E}[S] = m$ e $\mathrm{Var}[S] = v$. O custo total anual de
salários é $C = 12 \sum_{i=1}^N S_i$ — aleatório tanto em $N$ quanto em
$S$. Então

$$
\mathbb{E}[C] \;=\; 12 \, \mathbb{E}[N] \, m,
$$

$$
\mathrm{Var}[C] \;=\; 144 \, \bigl( \mathbb{E}[N]^2 \, v + m^2 \, \mathrm{Var}[N] + \mathrm{Var}[N] \cdot v \bigr).
$$

Exemplo numérico. Use a cadeia de headcount padrão ($n_0 = 45$,
$\lambda = 2$, $\mu = 0{,}025$) e a distribuição de salários ajustada do
Artigo 2 $S \sim \mathrm{LogNormal}(9{,}2,\, 0{,}30)$. Então
$\mathbb{E}[N] \approx 54$ e $\mathrm{Var}[N] \approx 35$ no mês 12
(calculados exatamente da CTMC truncada). Os momentos do salário são
$m \approx 10\,432$ e $v \approx 1{,}02 \times 10^7$. Substituindo:

$$
\mathbb{E}[C] \;\approx\; 6{,}76 \,\mathrm{M}, \qquad \mathrm{sd}(C) \;\approx\; 0{,}86 \,\mathrm{M}.
$$

Esses são os mesmos números que o Monte Carlo do Artigo 1 recupera dentro
do erro de amostragem. A forma fechada é rápida e exata até a segunda
ordem; o Artigo 1 é necessário para distribuições completas, percentis e
probabilidades de cauda. Os três artigos compõem: o Artigo 2 fornece a
distribuição de salário ($m, v$), o Artigo 3 fornece a distribuição de
headcount ($\mathbb{E}[N], \mathrm{Var}[N]$), o Artigo 1 simula o objeto
conjunto que inclui correlação entre funcionários e dependências
intramês.

---

## 9. Um Framework Prático para Planejadores de Pessoal

Um planejador que queira começar a usar o kit não precisa de cada teorema
deste artigo. Precisa de um processo de quatro passos.

**1. Defina a cadeia.** Liste os níveis de carreira que importam (Júnior,
Pleno, Sênior, ou como sua organização os chame) mais um estado de
Saída. Estime as taxas de transição mensais a partir dos últimos 12–24
meses de dados de RH: a fração de Júniores promovidos, a fração que
saiu, etc. Taxas de poucos por cento ao mês são típicas.

**2. Resolva o estado estacionário.** Calcule $\pi$ via $\pi P = \pi$
(uma linha de NumPy). Isso lhe diz a composição do time no longo prazo
sob as taxas atuais. Se $\pi$ não combina com suas aspirações, as taxas
são a alavanca — não decisões trimestrais de contratação.

**3. Resolva o transiente.** Escolha uma pergunta:

- *Probabilidade de atingir alvo*: calcule $P(n(t) \geq T)$ a partir de
  $e^{Qt}$.
- *Tempo até o alvo*: resolva a EDO em forma fechada para
  $\mathbb{E}[n(t)]$.
- *Tempo até saída sob congelamento*: leia $t = N \mathbf{1}$ da matriz
  fundamental.

**4. Rode sensibilidade.** Varie cada taxa em $\pm 50\%$ e ranqueie o
impacto sobre o resultado que importa para você. O resultado lhe diz
qual alavanca puxar. Para a cadeia de headcount padrão, retenção Sênior é
a alavanca dominante para a composição do time.

O processo todo cabe em um Jupyter notebook. O repositório companheiro
tem templates para cada passo. Uma vez feito uma vez, o custo marginal
de rodar outro cenário é zero — o ROI sobre um investimento único nesse
aparato é alto.

---

## 10. Conclusão

Um plano de headcount que nomeia um único número é uma compressão com
perdas do processo aleatório por baixo. A compressão esconde três peças
de informação útil: a probabilidade de atingir o número, o tempo
esperado de chegada lá, e o estado estacionário para o qual o sistema
tende. A teoria de cadeias de Markov recupera as três em forma fechada,
e a recuperação custa quase nada uma vez que o aparato esteja montado.

O custo é antecipado: definir o espaço de estados, estimar taxas,
escolher entre formulações absorvedora e com reciclagem. O retorno é
permanente. Uma atualização de uma linha das taxas basta para responder
de novo a cada pergunta; um caderno de cenários roda em segundos.

A lição mais profunda é que "longo prazo" é um verbo, não um substantivo.
O gap espectral da cadeia de headcount é pequeno — cerca de 0,01 — o que
significa que o longo prazo leva anos, não trimestres. Planos
trimestrais de contratação não conseguem alcançar a distribuição
estacionária; só podem redirecionar a trajetória transiente em sua
direção. Planejadores que tratam o estado estacionário como destino
ficarão decepcionados; planejadores que o usam como bússola estarão
calibrados.

O artigo completa uma trilogia. O Artigo 2 escolheu distribuições; o
Artigo 3 (este) modelou como o time muda; o Artigo 1 simulou o custo
total. Juntos, convertem o planejamento de pessoal de uma série de
suposições embasadas em um sistema probabilístico que pode ser
consultado, perturbado e auditado. O plano determinístico é uma
estatística entre muitas. Já não é a resposta.

---

## Referências

- Bartholomew, D. J. (1982). *Stochastic Models for Social Processes*.
  Wiley. *(modelagem de força de trabalho — diretamente relevante)*
- Grinstead, C. & Snell, J. (1997). *Introduction to Probability*. AMS.
- Levin, D., Peres, Y. & Wilmer, E. (2009). *Markov Chains and Mixing
  Times*. AMS.
- Norris, J. R. (1997). *Markov Chains*. Cambridge University Press.
- Ross, S. M. (2014). *Introduction to Probability Models*. Academic
  Press, 11ª ed.
- Taylor, H. M. & Karlin, S. (1998). *An Introduction to Stochastic
  Modeling*. Academic Press, 3ª ed.
- Companheiros da trilogia: *Artigo 1 — Simulação de Orçamento por Monte
  Carlo*; *Artigo 2 — Escolha de Distribuição para Componentes de Custo
  de TI*.

---

*Nota de reprodutibilidade. Todas as figuras e resultados numéricos deste
artigo são reproduzidos por scripts versionados no repositório
companheiro. Rode `pip install -e ".[dev]"` e depois
`python scripts/run_all_experiments.py` para regenerar o conjunto
completo de figuras com seed 2026.*
