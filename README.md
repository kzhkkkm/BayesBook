# ベイズ分析の理論と応用: R言語による経済データの分析

---

- [Rコード](#Rコード)
	- [データ](#データ)
	- [第2章](#第2章)
	- [第3章](#第3章)
	- [第4章](#第4章)
	- [第5章](#第5章)
	- [第6章](#第6章)
	- [付録A](#付録A)
	- [付録C](#付録C)
- [修正](#修正)
- [正誤表](#正誤表)

# Rコード

収束判定のための関数群: [MCMC.r](MCMC.r)

## データ
+ 為替レート: [ExRate.csv](Data/ExRate.csv)
+ 家計調査: [Income.csv](Data/Income.csv)
+ 工業統計: [mi.csv](Data/mi.csv)
+ 景気循環: [bc.csv](Data/bc.csv)

## 第2章
+ ヒストリカル・ボラティリティのべイズ分析: [Chap2-2.r](Chap2/Chap2-2.r)

## 第3章
+ 第3.1節 モンテカルロ法: [chap3-1.r](Chap3/Chap3-1.r)
+ 第3.3節 マルコフ連鎖モンテカルロ法: [Chap3-3.r](Chap3/Chap3-3.r)
+ 第3.3.1節 メトロポリス・ヘイスティングス・アルゴリズム: [Chap3-3-1.r](Chap3/Chap3-3-1.r)
+ 第3.3.2節 ギブズ・サンプラー: [Chap3-3-2.r](Chap3/Chap3-3-2.r)
+ 第3.3.3節 データ拡大法: [chap3-3-3.r](Chap3/Chap3-3-3.r)

## 第4章
+ 第4.1節 炭鉱爆発事故のベイズ分析: [Chap4-1.r](Chap4/Chap4-1.r)
+ 第4.2節 日本の所得分配の不平等のベイズ分析: [Chap4-2.r](Chap4/Chap4-2.r)
+ 第4.3節 円・ドルレートの為替収益率のベイズ分析: [Chap4-3.r](Chap4/Chap4-3.r)

## 第5章
+ 第5.1節 製造業における規模の経済性のベイズ分析: [Chap5-1.r](Chap5/Chap5-1.r)
+ 第5.2節 製造業における技術的効率性のベイズ分析: [Chap5-2.r](Chap5/Chap5-2.r)
+ 第5.3節 製造業における直接効果と間接効果のベイズ分析: [Chap5-3.r](Chap5/Chap5-3.r)

## 第6章
+ 第6.1節 景気循環のベイズ分析: [Chap6-1.r](Chap6/Chap6-1.r)

## 付録A
+ 付録A.2.7 最尤法:  [Chap-app-a.r](ChapAppendix/Chap-app-a.r)

## 付録C
+ 付録C.2 MHアルゴリズムの比較: [Chap-app-c.r](ChapAppendix/Chap-app-c.r)

---

# 修正
+ 図5.1の凡例を修正(22/11/16)
+ [chap3-3-3.r](Chap3/Chap3-3-3.r)75行目のバグを修正(22/11/17)
+ [chap3-3-1.r](Chap3/Chap3-3-1.r)，[chap3-3-3.r](Chap3/Chap3-3-3.r)，[chap4-2.r](Chap4/Chap4-2.r)，[Chap-app-c.r](ChapAppendix/Chap-app-c.r)を高速化のため修正(23/07/17)
---
# 正誤表
+ [正誤表](https://kzhkkkm.github.io/bayes.html)
