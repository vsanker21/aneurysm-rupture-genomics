"""
Analyze Supplementary PDFs
Extract information from downloaded supplementary PDFs
"""
import sys
from pathlib import Path

try:
    import PyPDF2
    PDF2_AVAILABLE = True
except ImportError:
    PDF2_AVAILABLE = False
    print("⚠️ PyPDF2 not available, trying alternative methods...")

try:
    import pdfplumber
    PDFPLUMBER_AVAILABLE = True
except ImportError:
    PDFPLUMBER_AVAILABLE = False

def analyze_pdf_pypdf2(pdf_path):
    """Analyze PDF using PyPDF2"""
    try:
        with open(pdf_path, 'rb') as pdf_file:
            reader = PyPDF2.PdfReader(pdf_file)
            
            info = {
                'pages': len(reader.pages),
                'title': reader.metadata.get('/Title', 'N/A') if reader.metadata else 'N/A',
                'subject': reader.metadata.get('/Subject', 'N/A') if reader.metadata else 'N/A',
                'author': reader.metadata.get('/Author', 'N/A') if reader.metadata else 'N/A',
                'first_page_text': ''
            }
            
            # Extract text from first page
            if len(reader.pages) > 0:
                try:
                    info['first_page_text'] = reader.pages[0].extract_text()[:1000]
                except:
                    info['first_page_text'] = 'Could not extract text'
            
            return info
    except Exception as e:
        return {'error': str(e)}

def analyze_pdf_pdfplumber(pdf_path):
    """Analyze PDF using pdfplumber"""
    try:
        with pdfplumber.open(pdf_path) as pdf:
            info = {
                'pages': len(pdf.pages),
                'first_page_text': ''
            }
            
            # Extract text from first page
            if len(pdf.pages) > 0:
                try:
                    info['first_page_text'] = pdf.pages[0].extract_text()[:1000] if pdf.pages[0].extract_text() else 'No text found'
                except:
                    info['first_page_text'] = 'Could not extract text'
            
            return info
    except Exception as e:
        return {'error': str(e)}

def analyze_pdfs():
    """Analyze all supplementary PDFs"""
    pdf_dir = Path("data/supplementary_pdfs")
    
    if not pdf_dir.exists():
        print(f"❌ Directory not found: {pdf_dir}")
        return
    
    pdfs = list(pdf_dir.glob("*.pdf"))
    
    if not pdfs:
        print(f"❌ No PDFs found in {pdf_dir}")
        return
    
    print("=" * 80)
    print("ANALYZING SUPPLEMENTARY PDFs")
    print("=" * 80)
    
    for pdf_path in pdfs:
        print(f"\n📄 Analyzing: {pdf_path.name}")
        print("-" * 80)
        
        # Try PyPDF2 first
        if PDF2_AVAILABLE:
            info = analyze_pdf_pypdf2(pdf_path)
            if 'error' not in info:
                print(f"  Pages: {info['pages']}")
                print(f"  Title: {info['title']}")
                print(f"  Subject: {info['subject']}")
                print(f"  Author: {info['author']}")
                print(f"\n  First page preview:")
                print(f"  {info['first_page_text'][:500]}...")
                continue
        
        # Try pdfplumber
        if PDFPLUMBER_AVAILABLE:
            info = analyze_pdf_pdfplumber(pdf_path)
            if 'error' not in info:
                print(f"  Pages: {info['pages']}")
                print(f"\n  First page preview:")
                print(f"  {info['first_page_text'][:500]}...")
                continue
        
        # Fallback: just file info
        print(f"  Size: {pdf_path.stat().st_size / (1024*1024):.2f} MB")
        print(f"  ⚠️ Could not extract text (install PyPDF2 or pdfplumber)")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    analyze_pdfs()






