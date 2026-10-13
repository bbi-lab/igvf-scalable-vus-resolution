#target illustrator
function cleanExport300() {
    if (app.documents.length === 0) {
        alert("Please open a document first!");
        return;
    }
    
    var doc = app.activeDocument;
    
    // 1. Create a dialog window to pick Format
    var dialog = new Window("dialog", "Select Clean Export Format");
    dialog.orientation = "column";
    dialog.alignChildren = ["left", "top"];
    
    dialog.add("statictext", undefined, "Choose your export format (Zero Prefixes):");
    var optPNG = dialog.add("radiobutton", undefined, "PNG (300 ppi / High Res)");
    var optPDF = dialog.add("radiobutton", undefined, "PDF (Individual Artboards)");
    optPNG.value = true; // default setting

    // PDF standard sub-choice (only relevant when optPDF is selected)
    var pdfStandardPanel = dialog.add("panel", undefined, "PDF Standard");
    pdfStandardPanel.orientation = "column";
    pdfStandardPanel.alignChildren = ["left", "top"];
    pdfStandardPanel.margins = 12;
    var optPDFStandard = pdfStandardPanel.add("radiobutton", undefined, "Standard (Acrobat 5)");
    var optPDFX3_2002 = pdfStandardPanel.add("radiobutton", undefined, "PDF/X-3:2002 (fixes missing/invisible placed PDFs)");
    var optPDFX3_2003 = pdfStandardPanel.add("radiobutton", undefined, "PDF/X-3:2003 (fixes missing/invisible placed PDFs)");
    optPDFX3_2003.value = true; // default to PDF/X-3 since it resolves the placed-PDF visibility bug

    var updatePdfStandardEnabled = function () {
        pdfStandardPanel.enabled = optPDF.value;
    };
    optPNG.onClick = updatePdfStandardEnabled;
    optPDF.onClick = updatePdfStandardEnabled;
    updatePdfStandardEnabled();

    var btnGroup = dialog.add("group");
    btnGroup.alignment = "right";
    var btnCancel = btnGroup.add("button", undefined, "Cancel", {name: "cancel"});
    var btnOK = btnGroup.add("button", undefined, "Export", {name: "ok"});
    
    if (dialog.show() == false) return; // Exit if user cancels
    
    // 2. Select the destination folder
    var targetFolder = Folder.selectDialog("Where do you want to save your clean files?");
    if (targetFolder == null) return;

    // 3. Process the export loop
    for (var i = 0; i < doc.artboards.length; i++) {
        doc.artboards.setActiveArtboardIndex(i);
        var abName = doc.artboards[i].name;
        
        if (optPNG.value) {
            // PNG Configured for 300 PPI via Scaling (416.66% scaling forces 72ppi up to 300ppi)
            var pngOptions = new ExportOptionsPNG24();
            pngOptions.antiAliasing = true;
            pngOptions.transparency = true;
            pngOptions.artBoardClipping = true;
            pngOptions.horizontalScale = 416.666; 
            pngOptions.verticalScale = 416.666;
            
            var pngFile = new File(targetFolder + "/" + abName + ".png");
            doc.exportFile(pngFile, ExportType.PNG24, pngOptions);
            
        } else if (optPDF.value) {
            // PDF Configuration
            var pdfOptions = new PDFSaveOptions();
            pdfOptions.generateAllBlurPoints = true;
            pdfOptions.artboardRange = (i + 1).toString(); // Target strictly this single artboard

            if (optPDFX3_2002.value || optPDFX3_2003.value) {
                // PDF/X-3 requires PDF 1.3 (Acrobat 4) compatibility -- matches Illustrator's
                // manual Save As dialog, which locks this field when a PDF/X-3 preset is chosen.
                // Leaving it at Acrobat 5 is what caused placed PDFs to still come out invisible.
                pdfOptions.compatibility = PDFCompatibility.ACROBAT4;
                pdfOptions.pDFXStandard = optPDFX3_2002.value ? PDFXStandard.PDFX32002 : PDFXStandard.PDFX32003;
                pdfOptions.preserveEditability = false;
            } else {
                pdfOptions.compatibility = PDFCompatibility.ACROBAT5;
            }

            var pdfFile = new File(targetFolder + "/" + abName + ".pdf");
            doc.saveAs(pdfFile, pdfOptions);
        }
    }
    
    alert("Export Complete! Checked and generated without document prefixes.");
}
cleanExport300();