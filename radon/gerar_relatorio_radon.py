import os
from radon.visitors import ComplexityVisitor
from radon.complexity import cc_visit
from radon.complexity import cc_rank
from radon.raw import analyze
from pathlib import Path
import html

# Caminhos
bioinf_path = Path(r"C:\Users\mende\OneDrive\Ambiente de Trabalho\U.M\Universidade\AASB\bioinf")
radon_path = Path(r"C:\Users\mende\OneDrive\Ambiente de Trabalho\U.M\Universidade\AASB\radon")
radon_path.mkdir(exist_ok=True)

# Arquivo HTML
html_file = radon_path / "relatorio_radon.html"

with open(html_file, "w", encoding="utf-8") as f:
    f.write("<html><head><meta charset='utf-8'><title>Relatório Radon</title></head><body>")
    f.write("<h1>Relatório Radon - Complexidade Ciclomática</h1>")

    # Percorre todos os arquivos Python
    for py_file in bioinf_path.rglob("*.py"):
        f.write(f"<h2>{py_file.name}</h2>")
        with open(py_file, "r", encoding="utf-8") as code_file:
            code = code_file.read()
            results = cc_visit(code)

            if not results:
                f.write("<p>Nenhuma função ou classe encontrada.</p>")
            else:
                f.write("<table border='1' cellpadding='5'><tr><th>Função/Classe</th><th>Complexidade</th><th>Rank</th></tr>")
                for item in results:
                    name = html.escape(item.name)
                    f.write(f"<tr><td>{name}</td><td>{item.complexity}</td><td>{cc_rank(item.complexity)}</td></tr>")
                f.write("</table>")

    f.write("</body></html>")

print(f"Relatório gerado em: {html_file}")
