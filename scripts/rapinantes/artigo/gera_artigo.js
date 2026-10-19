// gera_artigo.js
// Monta o artigo ABNT (.docx) a partir de reports/artigo/build/artigo_resolvido.md,
// que o referencias.py produz com as citações resolvidas e as referências em ::ref.
//
// Normas seguidas: NBR 6022:2018 (estrutura do artigo), NBR 6024:2012 (numeração
// progressiva), NBR 6028:2021 (resumo e palavras-chave), NBR 10520:2023 (citações)
// e NBR 6023:2018 (referências); geometria de página e tamanhos de fonte pela
// NBR 14724, que os modelos de artigo das universidades adotam; tabelas no padrão
// tabular do IBGE (sem fios verticais).
//
// Uso: node scripts/rapinantes/artigo/gera_artigo.js [entrada.md] [saida.docx]

"use strict";

const fs = require("fs");
const path = require("path");
const JSZip = require("jszip");
const {
  AlignmentType,
  BorderStyle,
  Document,
  Footer,
  FootnoteReferenceRun,
  Header,
  ImageRun,
  Packer,
  PageNumber,
  Paragraph,
  Table,
  TableCell,
  TableLayoutType,
  TableRow,
  TextRun,
  VerticalAlign,
  WidthType,
} = require("docx");

const RAIZ = path.resolve(__dirname, "..", "..", "..");
const ENTRADA = path.resolve(RAIZ, process.argv[2] || "reports/artigo/build/artigo_resolvido.md");
const SAIDA = path.resolve(
  RAIZ,
  process.argv[3] || "reports/artigo/Posicao_filogenetica_aves_de_rapina_ABNT.docx",
);

// 1 cm = 567 DXA. A4 com margens de 3 cm (superior e esquerda) e 2 cm (inferior e direita)
// deixa 16 cm de largura útil.
const CM = 567;
const LARGURA_UTIL = 16 * CM;
const FONTE = "Times New Roman";
const T12 = 24; // tamanhos em meio-ponto
const T10 = 20;
const LINHA_15 = 360; // entrelinha 1,5
const LINHA_1 = 240;
const RECUO = 709; // 1,25 cm
const IDIOMA = { value: "pt-BR" };

// ---------------------------------------------------------------- texto em linha

// **negrito**, *itálico* e ^sobrescrito^ podem se aninhar (ex.: [*S. l.*] dentro de
// uma referência); o marcador alterna o estado, e o texto entre marcadores herda
// os três estados vigentes.
// Espaço inseparável onde a quebra separaria o número da unidade ou do operador
// ("98,3 / %", "p = / 0,019", "< / 0,001").
function inseparaveis(texto) {
  return texto
    .replace(/(\d) %/g, "$1 %")
    .replace(/\bp = /g, "p = ")
    .replace(/([<>≥≤]) (\d)/g, "$1 $2");
}

function corridas(texto, base = {}) {
  texto = inseparaveis(texto);
  const runs = [];
  let negrito = false;
  let italico = false;
  let sobrescrito = false;
  let buffer = "";
  const descarrega = () => {
    if (!buffer) return;
    runs.push(
      new TextRun({
        text: buffer,
        font: FONTE,
        size: base.size || T12,
        bold: negrito || Boolean(base.bold),
        italics: italico !== Boolean(base.italics),
        superScript: sobrescrito,
        language: IDIOMA,
      }),
    );
    buffer = "";
  };
  for (let i = 0; i < texto.length; i++) {
    if (texto.startsWith("**", i)) {
      descarrega();
      negrito = !negrito;
      i++;
    } else if (texto[i] === "*") {
      descarrega();
      italico = !italico;
    } else if (texto[i] === "^") {
      descarrega();
      sobrescrito = !sobrescrito;
    } else {
      buffer += texto[i];
    }
  }
  descarrega();
  if (negrito || italico || sobrescrito) {
    throw new Error(`marcador sem par em: ${texto.slice(0, 80)}`);
  }
  return runs;
}

function semMarcadores(texto) {
  return texto.replace(/\*\*|\*|\^/g, "");
}

// ---------------------------------------------------------------- parágrafos

function corpo(texto) {
  return new Paragraph({
    children: corridas(texto),
    alignment: AlignmentType.JUSTIFIED,
    indent: { firstLine: RECUO },
    spacing: { line: LINHA_15, before: 0, after: 0 },
  });
}

function resumo(texto) {
  return new Paragraph({
    children: corridas(texto),
    alignment: AlignmentType.JUSTIFIED,
    spacing: { line: LINHA_1, before: 0, after: 0 },
  });
}

function palavrasChave(texto) {
  const [rotulo, ...resto] = texto.split(":");
  return new Paragraph({
    children: [
      new TextRun({ text: `${rotulo}:`, font: FONTE, size: T12, bold: true, language: IDIOMA }),
      ...corridas(resto.join(":")),
    ],
    alignment: AlignmentType.LEFT,
    spacing: { line: LINHA_1, before: LINHA_1, after: 0 },
  });
}

function titulo(nivel, texto) {
  // Seção primária em negrito e caixa alta, secundária em negrito, terciária em
  // itálico; os não numerados (resumo, referências) vão centralizados (NBR 6024).
  const opcoes = {
    keepNext: true,
    keepLines: true,
    spacing: { line: LINHA_15, before: LINHA_15, after: LINHA_15 },
    alignment: nivel === 0 ? AlignmentType.CENTER : AlignmentType.LEFT,
  };
  const estilo =
    nivel <= 1 ? { bold: true } : nivel === 2 ? { bold: true } : { italics: true };
  return new Paragraph({ ...opcoes, children: corridas(texto, estilo) });
}

function legenda(texto) {
  return new Paragraph({
    children: corridas(texto, { size: T10 }),
    alignment: AlignmentType.LEFT,
    keepNext: true,
    keepLines: true,
    spacing: { line: LINHA_1, before: LINHA_15, after: 80 },
  });
}

function fonteOuNota(texto, ultima) {
  return new Paragraph({
    children: corridas(texto, { size: T10 }),
    alignment: AlignmentType.JUSTIFIED,
    keepNext: !ultima,
    keepLines: true,
    spacing: { line: LINHA_1, before: 60, after: ultima ? LINHA_15 : 0 },
  });
}

function referencia(texto) {
  // NBR 6023: espaço simples, alinhada à esquerda, separadas por uma linha em branco.
  return new Paragraph({
    children: corridas(texto),
    alignment: AlignmentType.LEFT,
    spacing: { line: LINHA_1, before: 0, after: LINHA_1 },
  });
}

// ---------------------------------------------------------------- tabela

// Cada coluna recebe ao menos a largura da sua palavra mais longa (senão o Word e o
// LibreOffice quebram "Veredito" no meio); a sobra vai para as colunas na proporção
// do volume de texto, sem passar do que a coluna precisa para caber numa linha.
const DXA_POR_CARACTERE = 105; // média do Times New Roman 10 pt, com folga para o negrito
const MARGEM_CELULA = 160;

function larguras(linhas) {
  const ncol = linhas[0].length;
  const minimo = [];
  const desejado = [];
  const volume = [];
  for (let c = 0; c < ncol; c++) {
    const textos = linhas.map((l) => semMarcadores(l[c] || ""));
    const palavra = Math.max(...textos.flatMap((t) => t.split(/\s+/).map((p) => p.length)));
    const linha = Math.max(...textos.map((t) => t.length));
    minimo.push(palavra * DXA_POR_CARACTERE + MARGEM_CELULA);
    desejado.push(Math.max(linha * DXA_POR_CARACTERE + MARGEM_CELULA, minimo[c]));
    volume.push(textos.reduce((a, t) => a + t.length, 0));
  }
  const soma = (v) => v.reduce((a, b) => a + b, 0);
  let cols;
  if (soma(desejado) <= LARGURA_UTIL) {
    cols = desejado.map((d) => (d / soma(desejado)) * LARGURA_UTIL);
  } else if (soma(minimo) >= LARGURA_UTIL) {
    cols = minimo.map((m) => (m / soma(minimo)) * LARGURA_UTIL);
  } else {
    // Enchimento por volume: a coluna que atinge o desejado sai da partilha e a sobra
    // dela volta para as outras.
    cols = [...minimo];
    let abertas = cols.map((_, c) => c);
    while (abertas.length) {
      const sobra = LARGURA_UTIL - soma(cols);
      const pesoTotal = soma(abertas.map((c) => volume[c]));
      if (sobra < 1 || pesoTotal === 0) break;
      const cheias = [];
      for (const c of abertas) {
        cols[c] += (volume[c] / pesoTotal) * sobra;
        if (cols[c] >= desejado[c]) {
          cols[c] = desejado[c];
          cheias.push(c);
        }
      }
      if (!cheias.length) break;
      abertas = abertas.filter((c) => !cheias.includes(c));
    }
  }
  cols = cols.map(Math.floor);
  cols[cols.length - 1] += LARGURA_UTIL - cols.reduce((a, b) => a + b, 0);
  return cols;
}

function tabela(linhas) {
  const cols = larguras(linhas);
  const nada = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
  const fio = { style: BorderStyle.SINGLE, size: 6, color: "000000" };
  const ultima = linhas.length - 1;
  const rows = linhas.map(
    (celulas, r) =>
      new TableRow({
        tableHeader: r === 0,
        cantSplit: true,
        children: celulas.map(
          (texto, c) =>
            new TableCell({
              width: { size: cols[c], type: WidthType.DXA },
              verticalAlign: VerticalAlign.TOP,
              margins: { top: 40, bottom: 40, left: 80, right: 80 },
              // Padrão IBGE: fio no topo e na base do cabeçalho e na base da tabela.
              borders: {
                top: r === 0 ? fio : nada,
                bottom: r === 0 || r === ultima ? fio : nada,
                left: nada,
                right: nada,
              },
              children: [
                // keepNext em todas as linhas mantém a tabela inteira na mesma página
                // e presa à Fonte logo abaixo; nenhuma tabela do artigo passa de meia página.
                new Paragraph({
                  children: corridas(texto, { size: T10, bold: r === 0 }),
                  alignment: AlignmentType.LEFT,
                  keepNext: true,
                  keepLines: true,
                  spacing: { line: LINHA_1, before: 0, after: 0 },
                }),
              ],
            }),
        ),
      }),
  );
  return new Table({
    rows,
    width: { size: LARGURA_UTIL, type: WidthType.DXA },
    columnWidths: cols,
    layout: TableLayoutType.FIXED,
    borders: {
      top: nada,
      bottom: nada,
      left: nada,
      right: nada,
      insideHorizontal: nada,
      insideVertical: nada,
    },
  });
}

function celulasDe(linha) {
  return linha
    .trim()
    .replace(/^\|/, "")
    .replace(/\|$/, "")
    .split("|")
    .map((c) => c.trim());
}

// ---------------------------------------------------------------- figura

function dimensoesPng(buf) {
  if (buf.toString("ascii", 1, 4) !== "PNG") throw new Error("figura não é PNG");
  return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) };
}

function figura(caminho, larguraCm) {
  const buf = fs.readFileSync(path.resolve(RAIZ, caminho));
  const { w, h } = dimensoesPng(buf);
  const px = (larguraCm / 2.54) * 96;
  const altura = (px * h) / w;
  const alturaCm = (altura / 96) * 2.54;
  if (alturaCm > 22) throw new Error(`${caminho}: ${alturaCm.toFixed(1)} cm não cabe na página`);
  return new Paragraph({
    children: [
      new ImageRun({
        type: "png",
        data: buf,
        transformation: { width: Math.round(px), height: Math.round(altura) },
        altText: { title: path.basename(caminho), description: path.basename(caminho), name: path.basename(caminho) },
      }),
    ],
    alignment: AlignmentType.CENTER,
    keepNext: true,
    // Sem entrelinha explícita: com w:line o LibreOffice corta a imagem na altura
    // de uma linha de texto.
    spacing: { before: 0, after: 0 },
  });
}

// ---------------------------------------------------------------- leitura

function lerFrontMatter(texto) {
  const m = texto.match(/^---\n([\s\S]*?)\n---\n/);
  if (!m) throw new Error("front matter ausente");
  const meta = {};
  for (const linha of m[1].split("\n")) {
    const i = linha.indexOf(":");
    meta[linha.slice(0, i).trim()] = linha.slice(i + 1).trim();
  }
  return { meta, resto: texto.slice(m[0].length) };
}

function contarPalavras(texto) {
  return semMarcadores(texto).split(/\s+/).filter((p) => /[\p{L}\p{N}]/u.test(p)).length;
}

function montar(texto) {
  const { meta, resto } = lerFrontMatter(texto.replace(/\r\n/g, "\n"));
  const filhos = [];

  // Página de abertura: título: subtítulo, título em inglês, autoria com nota de rodapé.
  filhos.push(
    new Paragraph({
      children: [
        ...corridas(meta.titulo, { bold: true }),
        ...corridas(`: ${meta.subtitulo}`),
      ],
      alignment: AlignmentType.CENTER,
      spacing: { line: LINHA_15, before: 0, after: 0 },
    }),
    new Paragraph({
      children: corridas(`${meta.titulo_en}: ${meta.subtitulo_en}`, { italics: true }),
      alignment: AlignmentType.CENTER,
      spacing: { line: LINHA_15, before: LINHA_15, after: 0 },
    }),
    new Paragraph({
      children: [...corridas(meta.autor), new FootnoteReferenceRun(1)],
      alignment: AlignmentType.RIGHT,
      spacing: { line: LINHA_15, before: LINHA_15 * 2, after: 0 },
    }),
  );

  const linhas = resto.split("\n");
  let secao = "";
  let tabelaAberta = null;
  const avisos = [];

  const fecharTabela = () => {
    if (!tabelaAberta) return;
    const semSeparador = tabelaAberta.filter((l) => !/^\|?\s*:?-{3,}/.test(l.trim()));
    filhos.push(tabela(semSeparador.map(celulasDe)));
    tabelaAberta = null;
  };

  for (let i = 0; i < linhas.length; i++) {
    const linha = linhas[i];
    const proxima = (linhas[i + 1] || "").trim();

    if (tabelaAberta && linha.trim().startsWith("|")) {
      tabelaAberta.push(linha);
      continue;
    }
    fecharTabela();

    if (!linha.trim()) continue;

    if (linha.startsWith("#* ")) {
      secao = linha.slice(3).trim();
      filhos.push(titulo(0, secao));
    } else if (/^#{1,3} /.test(linha)) {
      const nivel = linha.match(/^#+/)[0].length;
      secao = linha.replace(/^#+ /, "").trim();
      filhos.push(titulo(nivel, secao));
    } else if (linha.startsWith("::tabela ")) {
      filhos.push(legenda(linha.slice(9).trim()));
      tabelaAberta = [];
    } else if (linha.startsWith("::figura ")) {
      const [caminho, largura, texto] = linha.slice(9).split("|");
      filhos.push(legenda(texto.trim()), figura(caminho.trim(), Number(largura)));
    } else if (linha.startsWith("::fonte ") || linha.startsWith("::nota ")) {
      const ultima = !proxima.startsWith("::nota ");
      filhos.push(fonteOuNota(linha.replace(/^::(fonte|nota) /, ""), ultima));
    } else if (linha.startsWith("::ref ")) {
      filhos.push(referencia(linha.slice(6).trim()));
    } else if (linha.startsWith("::")) {
      throw new Error(`marcador desconhecido: ${linha.slice(0, 40)}`);
    } else if (/^(Palavras-chave|Keywords):/.test(linha)) {
      filhos.push(palavrasChave(linha.trim()));
      const termos = linha.split(":")[1].trim();
      const n = termos.replace(/\.$/, "").split(";").length;
      if (!termos.endsWith(".") || n < 3 || n > 5) avisos.push(`${linha.split(":")[0]}: ${n} termos`);
      if (linha.startsWith("Keywords")) {
        // NBR 6022:2018, 5.1.6: datas de submissão e aprovação, preenchidas pelo periódico.
        filhos.push(
          new Paragraph({
            children: corridas("Artigo submetido em: [a preencher]. Aprovado em: [a preencher]."),
            alignment: AlignmentType.LEFT,
            spacing: { line: LINHA_1, before: LINHA_15, after: 0 },
          }),
        );
      }
    } else if (secao === "RESUMO" || secao === "ABSTRACT") {
      const n = contarPalavras(linha);
      if (n < 100 || n > 250) avisos.push(`${secao}: ${n} palavras (NBR 6028 pede 100 a 250)`);
      console.log(`${secao}: ${n} palavras`);
      filhos.push(resumo(linha.trim()));
    } else {
      filhos.push(corpo(linha.trim()));
    }
  }
  fecharTabela();
  return { meta, filhos, avisos };
}

// ---------------------------------------------------------------- documento

async function principal() {
  const { meta, filhos, avisos } = montar(fs.readFileSync(ENTRADA, "utf8"));

  const doc = new Document({
    creator: meta.autor,
    title: `${semMarcadores(meta.titulo)}: ${meta.subtitulo}`,
    description: "Artigo científico nas normas da ABNT (NBR 6022:2018).",
    language: "pt-BR",
    styles: {
      default: {
        document: { run: { font: FONTE, size: T12, language: IDIOMA } },
      },
    },
    footnotes: {
      1: {
        children: [
          new Paragraph({
            children: corridas(` ${meta.nota_autor}`, { size: T10 }),
            alignment: AlignmentType.JUSTIFIED,
            spacing: { line: LINHA_1, before: 0, after: 0 },
          }),
        ],
      },
    },
    sections: [
      {
        properties: {
          page: {
            size: { width: 11906, height: 16838 },
            margin: { top: 3 * CM, left: 3 * CM, bottom: 2 * CM, right: 2 * CM, header: 2 * CM, footer: CM },
          },
        },
        headers: {
          default: new Header({
            children: [
              new Paragraph({
                alignment: AlignmentType.RIGHT,
                children: [new TextRun({ children: [PageNumber.CURRENT], font: FONTE, size: T10 })],
              }),
            ],
          }),
        },
        footers: { default: new Footer({ children: [new Paragraph({ children: [] })] }) },
        children: filhos,
      },
    ],
  });

  // A NBR 6022 indica a nota de autoria por asterisco; o docx-js só numera notas,
  // então a marca personalizada entra direto no XML.
  const zip = await JSZip.loadAsync(await Packer.toBuffer(doc));
  const documento = await zip.file("word/document.xml").async("string");
  const notas = await zip.file("word/footnotes.xml").async("string");
  const ref = /<w:footnoteReference w:id="1"\/>/;
  if (!ref.test(documento) || !notas.includes("<w:footnoteRef/>")) {
    throw new Error("nota de autoria não encontrada no XML gerado");
  }
  zip.file(
    "word/document.xml",
    documento.replace(ref, '<w:footnoteReference w:customMarkFollows="1" w:id="1"/><w:t>*</w:t>'),
  );
  zip.file("word/footnotes.xml", notas.replace("<w:footnoteRef/>", "<w:t>*</w:t>"));
  fs.writeFileSync(SAIDA, await zip.generateAsync({ type: "nodebuffer", compression: "DEFLATE" }));

  console.log(path.relative(RAIZ, SAIDA));
  for (const a of avisos) console.log(`AVISO ${a}`);
  if (avisos.length) process.exitCode = 1;
}

principal().catch((e) => {
  console.error(e.message);
  process.exit(2);
});
