BRUNO B BIZARRIA
18/02/2021

CÓDIGO PARA OCUPAR HALOS COM GALÁXIAS SIMULANDO O DES

O código povoa um lightcone de halos com galáxias, utilizando uma função do número médio de galáxias dentro de um halo
dado a massa do halo. Esta função é obtida por observação, e à princípio seria através dela que simularíamos nosso survey.
É necessário considerar algum tipo de corte em magnitude, assim como o survey observa. Uma tentativa seria considerar
uma população de galáxias observadas delimitadas em redsfhifts diferentes, e obter os parâmetros de fit da função para 
cada uma dessas amostras. 

Vamos ao presente caso:

Estamos considerando o DES, com as amostras RedMagic, nos bins [0, 0.35] e [0.35, 0.50]. Então temos duas funções de povoamento, obtidas por  Zacharegkas(arxiv:2106.08438).
Portanto, como se trata de uma versão beta do código, é necessário rodar dois códigos de povoamento, para os dois intervalos 
de redshift. O erro é aplicado para cada galáxia, sendo que isto ocorre no STEP 4. Como produto, são gerados duas pastas, cada qual com uma pasta densidade e outra de contagem. A primeira é densidade numérica, e a outra é número absoluto, ambas com 30 canais.

procedimento para rodar os códigos:

1-salvar simulador.py e def_funcs.py na mesma pasta

2-Verificar se no STEP2 está o lightcone total ou o parcial

3-simulador.py Colocar como path o local onde você deseja salvar as pastas com os arquivos intermediários e os catálogos finais

4-Os catálogos finais não foram feitos, que é a superposição dos dois samples. Preferi deixar assim para facilitar na comparação com o DES

Catálogos de saída:
Os catálogos resultantes tem seu nome no formato nMhz.npy. A relação entre n e o intervalo de redshift que ele representa está explícita abaixo:
0	0.449-0.435
1	0.435-0.422
2	0.422-0.409
3	0.409-0.396
4	0.396-0.383
5	0.383-0.371
6	0.371-0.358
7	0.358-0.346
8	0.346-0.334
9	0.334-0.323
10	0.323-0.311
11	0.311-0.300
12	0.300-0.289
13	0.289-0.278
14	0.278-0.268
15	0.268-0.257
16	0.257-0.247
17	0.247-0.237
18	0.237-0.227
19	0.227-0.217
20	0.217-0.207
21	0.207-0.198
22	0.198-0.188
23	0.188-0.179
24	0.179-0.170
25	0.170-0.161
26	0.161-0.152
27	0.152-0.144
28	0.144-0.135
29	0.135-0.127

Dois tipos de catálogos são gerados. Um de número de galáxias por pixel (pasta contagem), e o outro de densidade de galáxias por pixel (densidade), que é o número de galáxias por Mpc^3 correspondente ao volume que o pixel representa.
