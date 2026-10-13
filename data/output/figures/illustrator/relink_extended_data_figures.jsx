#target illustrator

// Renumbers linked/placed files whose path contains an
// "extended_data_figure_N" directory segment, per the mapping below, and
// relinks the active document to the new paths. Suffixes on the directory
// name (e.g. "extended_data_figure_9alt", "extended_data_figure_6alt_missense",
// "extended_data_figure_9_op") are preserved; only the number is remapped.
// Artboards whose names contain "extended_data_figure_N" (or the
// human-readable "Extended Data Figure N") are renamed the same way.
//
// Covers both linked PlacedItem objects (placed PDFs/AI/EPS) and linked
// RasterItem objects (placed PNG/TIFF/etc. that were not embedded).
//
// Three complications this script works around:
//
// 1. Many PageItems in a typical assembled figure (drop shadows, rasterized
//    effects, clipping masks, etc.) report as RasterItem/PlacedItem but have
//    no backing file at all -- reading `.file` on one of those throws "There
//    is no file associated with this item". That's expected and not a sign
//    anything is wrong; those items are skipped quietly.
//
// 2. If the on-disk asset folders have ALREADY been renamed (e.g. you ran
//    the figure-renumbering refactor before relinking Illustrator), every
//    real link is now "missing" in the Links panel -- and for a missing
//    link, `.file` ALSO throws that same "no file associated" error, even
//    though the Links panel still shows its original filename. In that case
//    this script falls back to searching the new asset tree (assumed to
//    live at "<this .ai file's folder>/../assets", i.e.
//    data/output/figures/assets/ in this repo) for a file with the same
//    basename. Filenames themselves never embed the figure number, only
//    their containing directory does -- but the SAME basename can
//    legitimately exist under more than one figure's directory (e.g. a
//    "main" and "alt" version of the same panel, or two unrelated figures
//    that happen to reuse a generic name), so a basename match alone can be
//    ambiguous. When that happens, the script narrows the candidates using
//    the artboard the broken item sits on: each artboard corresponds to one
//    figure, so if exactly one candidate's renumbered figure identifier
//    (number + alt/missense/op-style suffix, parsed the same way artboard
//    names are) matches the artboard's own (pre-rename) identifier, that
//    candidate wins. If more than one candidate still matches (or none
//    does), it's reported as ambiguous rather than guessed.
//
// 3. Despite Adobe's own scripting reference describing PlacedItem.file as
//    read/write, reassigning `.file` on an EXISTING placed item throws
//    "Unable to set placed item's file" in current Illustrator -- it only
//    works on a freshly created one. So relinking is done by creating a new
//    PlacedItem in the same parent (layer/group), pointing it at the new
//    file, copying the old item's transform matrix so position/size/
//    rotation match exactly, inserting it at the old item's z-order slot,
//    and then removing the old item.

var RENUMBER_MAP = {
    4: 2,
    6: 3,
    7: 4,
    8: 5,
    9: 6,
    10: 7,
    11: 8,
    12: 9,
    16: 10
};

var DIR_PATTERN = /^extended_data_figure_(\d+)(.*)$/;
var ARTBOARD_PATTERN = /(extended[_ ]data[_ ]figure[_ ]?)(\d+)([A-Za-z_]*)/gi;
var IDENTIFIER_PATTERN = /extended[_ ]data[_ ]figure[_ ]?(\d+)[_ ]?(.*)$/i;

// "_alt", "_Alt_Missense", "alt_missense" all look the same once separators
// are stripped and case is flattened -- this is what lets an artboard name
// like "extended_data_figure_6_alt_missense" compare against a directory
// name like "extended_data_figure_3alt_missense".
function normalizeSuffix(s) {
    return s.replace(/[^a-zA-Z0-9]/g, "").toLowerCase();
}

// Parses "extended_data_figure_6_alt_missense_captioned" (or the directory-
// style "extended_data_figure_3alt_missense") into {number, suffix}, with
// suffix normalized. Returns null if the name doesn't match at all.
function parseFigureIdentifier(name) {
    var m = IDENTIFIER_PATTERN.exec(name);
    if (!m) return null;
    return { number: parseInt(m[1], 10), suffix: normalizeSuffix(m[2]) };
}

function mapPathSegment(segment) {
    var m = DIR_PATTERN.exec(segment);
    if (!m) return null;
    var oldNum = parseInt(m[1], 10);
    var suffix = m[2];
    if (!(oldNum in RENUMBER_MAP)) return null;
    return "extended_data_figure_" + RENUMBER_MAP[oldNum] + suffix;
}

function mapArtboardName(name) {
    var changed = false;
    var newName = name.replace(ARTBOARD_PATTERN, function (full, prefix, num, suffix) {
        var oldNum = parseInt(num, 10);
        if (!(oldNum in RENUMBER_MAP)) return full;
        changed = true;
        return prefix + RENUMBER_MAP[oldNum] + suffix;
    });
    return changed ? newName : null;
}

function remapPath(fsName) {
    var sep = (fsName.indexOf("\\") !== -1 && fsName.indexOf("/") === -1) ? "\\" : "/";
    var parts = fsName.split(sep);
    var changed = false;
    for (var i = 0; i < parts.length; i++) {
        var mapped = mapPathSegment(parts[i]);
        if (mapped !== null) {
            parts[i] = mapped;
            changed = true;
        }
    }
    if (!changed) return null;
    return parts.join(sep);
}

// Recursively indexes every file under `root` by basename, so a broken
// link's filename (still available even when its full path isn't) can be
// matched against the already-renumbered asset tree.
function buildFileIndex(root) {
    var index = {};
    function walk(folder) {
        var entries;
        try {
            entries = folder.getFiles();
        } catch (e) {
            return;
        }
        for (var i = 0; i < entries.length; i++) {
            var entry = entries[i];
            if (entry instanceof Folder) {
                walk(entry);
            } else if (entry instanceof File) {
                if (!index[entry.name]) index[entry.name] = [];
                index[entry.name].push(entry);
            }
        }
    }
    walk(root);
    return index;
}

// Relinking by reassigning `.file` on an existing PlacedItem fails in
// current Illustrator ("Unable to set placed item's file"); only setting
// `.file` right after creating a new PlacedItem works. So: create a new
// item in the same parent, point it at the new file, copy the old item's
// matrix (position/size/rotation) and name, and splice it into the old
// item's z-order slot before removing the old one.
function relinkByRecreate(oldItem, newPath) {
    var parent = oldItem.parent;
    var matrix = oldItem.matrix;
    var oldName = null;
    try { oldName = oldItem.name; } catch (eName) { oldName = null; }

    var newItem = parent.placedItems.add();
    try {
        newItem.file = new File(newPath);
        newItem.matrix = matrix;
        newItem.move(oldItem, ElementPlacement.PLACEBEFORE);
        if (oldName) {
            try { newItem.name = oldName; } catch (eRename) { /* non-fatal */ }
        }
    } catch (e) {
        try { newItem.remove(); } catch (eCleanup) { /* best effort */ }
        throw e;
    }
    oldItem.remove();
    return newItem;
}

// Captures each artboard's rect and parsed (pre-rename) figure identifier,
// for later use disambiguating broken links by which artboard they sit on.
// Must be called BEFORE any artboard renaming happens.
function captureArtboardInfo(doc) {
    var info = [];
    for (var i = 0; i < doc.artboards.length; i++) {
        var ab = doc.artboards[i];
        var identifier = parseFigureIdentifier(ab.name);
        if (identifier !== null) {
            info.push({ rect: ab.artboardRect, identifier: identifier, name: ab.name });
        }
    }
    return info;
}

// Finds the single artboard whose rect contains the item's center point, and
// returns its parsed (old) figure identifier, or null if the item's bounds
// aren't readable or it doesn't land uniquely on one of the known artboards.
function findArtboardIdentifierForItem(item, artboardInfo) {
    var bounds;
    try {
        bounds = item.geometricBounds; // [left, top, right, bottom]
    } catch (e) {
        return null;
    }
    var cx = (bounds[0] + bounds[2]) / 2;
    var cy = (bounds[1] + bounds[3]) / 2;
    var match = null;
    for (var i = 0; i < artboardInfo.length; i++) {
        var r = artboardInfo[i].rect; // [left, top, right, bottom]; top > bottom
        if (cx >= r[0] && cx <= r[2] && cy <= r[1] && cy >= r[3]) {
            if (match !== null) return null; // sits on more than one -- ambiguous itself
            match = artboardInfo[i];
        }
    }
    return match ? match.identifier : null;
}

// Given a broken link's basename candidates (files, not yet disambiguated)
// and the renumbered identifier an artboard implies, returns the single
// candidate matching that identifier, or null if that's still not unique.
function disambiguateByArtboard(candidates, artboardIdentifier) {
    if (artboardIdentifier === null) return null;
    if (!(artboardIdentifier.number in RENUMBER_MAP)) return null;
    var expectedNumber = RENUMBER_MAP[artboardIdentifier.number];
    var expectedSuffix = artboardIdentifier.suffix;

    var matches = [];
    for (var i = 0; i < candidates.length; i++) {
        var candIdentifier = parseFigureIdentifier(candidates[i].fsName);
        if (!candIdentifier || candIdentifier.number !== expectedNumber) continue;
        var candSuffix = candIdentifier.suffix;
        // A bare (no-suffix) candidate only matches a bare artboard identifier;
        // a suffixed candidate matches if its suffix is part of the artboard's
        // (which may carry extra decorative tokens like "captioned").
        var ok = candSuffix === "" ? expectedSuffix === "" : expectedSuffix.indexOf(candSuffix) !== -1;
        if (ok) matches.push(candidates[i]);
    }
    return matches.length === 1 ? matches[0] : null;
}

function collectLinkedItems(doc) {
    var items = [];
    var i;
    for (i = 0; i < doc.placedItems.length; i++) {
        items.push({ kind: "PlacedItem", item: doc.placedItems[i] });
    }
    for (i = 0; i < doc.rasterItems.length; i++) {
        var r = doc.rasterItems[i];
        var embedded = true;
        try {
            embedded = r.embedded;
        } catch (eEmbedded) {
            // Can't tell -- treat like anything else with no readable file.
        }
        if (!embedded) items.push({ kind: "RasterItem", item: r });
    }
    return items;
}

function relinkExtendedDataFigures() {
    if (app.documents.length === 0) {
        alert("Please open the Illustrator document first!");
        return;
    }

    var dialog = new Window("dialog", "Relink Extended Data Figure Assets");
    dialog.orientation = "column";
    dialog.alignChildren = ["left", "top"];
    dialog.add("statictext", undefined, "Renumber extended_data_figure_N asset folders in linked files");
    dialog.add("statictext", undefined, "and rename matching artboards.");
    var dryRunCheck = dialog.add("checkbox", undefined, "Dry run (preview only, make no changes)");
    dryRunCheck.value = true;

    var btnGroup = dialog.add("group");
    btnGroup.alignment = "right";
    btnGroup.add("button", undefined, "Cancel", {name: "cancel"});
    btnGroup.add("button", undefined, "Run", {name: "ok"});

    if (dialog.show() == false) return;
    var dryRun = dryRunCheck.value;

    var doc = app.activeDocument;
    var items = collectLinkedItems(doc);
    // Must happen before any artboard renaming below, since it needs each
    // artboard's ORIGINAL (pre-rename) name to parse its old figure identifier.
    var artboardInfo = captureArtboardInfo(doc);

    // Assumes this .ai file lives directly under data/output/figures/Illustrator/
    // (as every .ai in this repo does), so its asset tree is a sibling of
    // that folder: data/output/figures/assets/.
    var assetsRoot = new Folder(doc.fullName.parent.parent.fsName + "/assets");
    var fileIndex = assetsRoot.exists ? buildFileIndex(assetsRoot) : {};

    var changedCount = 0;
    var skipCount = 0;
    var noFileCount = 0;
    var failCount = 0;
    var ambiguousCount = 0;
    var notFoundCount = 0;
    var report = [];
    var fileLog = [];

    for (var a = 0; a < doc.artboards.length; a++) {
        var artboard = doc.artboards[a];
        var oldName = artboard.name;
        var newName = mapArtboardName(oldName);
        if (newName === null) {
            skipCount++;
            continue;
        }
        if (dryRun) {
            report.push("WOULD RENAME ARTBOARD:\n  \"" + oldName + "\"\n  -> \"" + newName + "\"");
            changedCount++;
            continue;
        }
        try {
            artboard.name = newName;
            report.push("RENAMED ARTBOARD:\n  \"" + oldName + "\"\n  -> \"" + newName + "\"");
            changedCount++;
        } catch (eArtboard) {
            failCount++;
            report.push("FAILED to rename artboard (" + eArtboard + "):\n  \"" + oldName + "\"\n  -> \"" + newName + "\"");
        }
    }

    for (var i = 0; i < items.length; i++) {
        var kind = items[i].kind;
        var item = items[i].item;
        var oldFile = null;
        try {
            oldFile = item.file;
        } catch (eRead) {
            oldFile = null;
        }

        var oldPath, newPath;

        if (oldFile !== null) {
            // Link resolves fine -- remap its full recorded path directly.
            oldPath = oldFile.fsName;
            fileLog.push(kind + " " + i + ": " + oldPath + (oldFile.exists ? "" : "  [FILE MISSING ON DISK]"));
            newPath = remapPath(oldPath);
            if (newPath === null) {
                skipCount++;
                continue;
            }
        } else {
            // Broken/missing link (or no file at all) -- .file throws either
            // way. Try the item's bare name and search the new asset tree.
            var baseName = null;
            try {
                baseName = item.name;
            } catch (eName) {
                baseName = null;
            }
            if (!baseName) {
                noFileCount++;
                continue;
            }
            var candidates = fileIndex[baseName];
            if (!candidates || candidates.length === 0) {
                notFoundCount++;
                fileLog.push(kind + " " + i + ": [BROKEN LINK] \"" + baseName + "\" -- no match found under " + assetsRoot.fsName);
                continue;
            }
            if (candidates.length > 1) {
                var artboardIdentifier = findArtboardIdentifierForItem(item, artboardInfo);
                var resolved = disambiguateByArtboard(candidates, artboardIdentifier);
                if (resolved === null) {
                    ambiguousCount++;
                    var candidateList = [];
                    for (var c = 0; c < candidates.length; c++) candidateList.push(candidates[c].fsName);
                    fileLog.push(kind + " " + i + ": [BROKEN LINK] \"" + baseName + "\" -- AMBIGUOUS, " + candidates.length + " matches" +
                        (artboardIdentifier ? " (artboard didn't narrow it down)" : " (couldn't determine the item's artboard)") +
                        ":\n    " + candidateList.join("\n    "));
                    continue;
                }
                oldPath = "[broken link] " + baseName;
                newPath = resolved.fsName;
                fileLog.push(kind + " " + i + ": [BROKEN LINK] \"" + baseName + "\" -- " + candidates.length + " basename matches, resolved via its artboard to " + newPath);
            } else {
                oldPath = "[broken link] " + baseName;
                newPath = candidates[0].fsName;
                fileLog.push(kind + " " + i + ": [BROKEN LINK] \"" + baseName + "\" -- resolved by name to " + newPath);
            }
        }

        if (dryRun) {
            report.push("WOULD RELINK:\n  " + oldPath + "\n  -> " + newPath);
            changedCount++;
            continue;
        }
        try {
            relinkByRecreate(item, newPath);
            report.push("RELINKED:\n  " + oldPath + "\n  -> " + newPath);
            changedCount++;
        } catch (eWrite) {
            failCount++;
            report.push("FAILED to relink (" + eWrite + "):\n  " + oldPath + "\n  -> " + newPath);
        }
    }

    var summary = (dryRun ? "[DRY RUN] " : "") +
        "Relink/rename summary for \"" + doc.name + "\"\n" +
        "  Matched & " + (dryRun ? "would change" : "changed") + " (links + artboards): " + changedCount + "\n" +
        "  Items with a resolvable file, but no mapped figure number: " + skipCount + "\n" +
        "  Broken links resolved by filename search: see log\n" +
        "  Broken links with NO matching filename in " + assetsRoot.fsName + ": " + notFoundCount + "\n" +
        "  Broken links with AMBIGUOUS (multiple) filename matches: " + ambiguousCount + "\n" +
        "  Items with no file at all (effects/rasterized art -- expected, not an error): " + noFileCount + "\n" +
        "  Failed: " + failCount;

    var logFile = new File(doc.fullName.fsName.replace(/\.ai$/i, "") + "_relink_log.txt");
    logFile.open("w");
    logFile.write(
        summary + "\n\n" +
        "Asset search root: " + assetsRoot.fsName + (assetsRoot.exists ? "" : "  [DOES NOT EXIST]") + "\n\n" +
        "--- Every linked item with a file or a broken-link name (for diagnosing matches) ---\n" +
        fileLog.join("\n") + "\n\n" +
        "--- Changes ---\n" +
        report.join("\n\n")
    );
    logFile.close();

    alert(summary + "\n\nFull log written to:\n" + logFile.fsName +
        (dryRun ? "\n\nUncheck \"Dry run\" and run again to apply." :
                  "\n\nReview the Links panel, then save the document."));
}

relinkExtendedDataFigures();
