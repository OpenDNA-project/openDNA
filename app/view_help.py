"""What the "?" and the cog next to a pane's view picker show, one entry per view of the picker (workspace.py: VIEWS and
TAB_TOOLS). Edit the texts here; tests/test_view_help.py checks that every view has one.

HELP[view]    - the short guide: 5-8 bullet lines, plain language (keep it near 120 words).
Every key gesture below is one the code handles (see also shortcuts.py, which lists them for Help > Keyboard Shortcuts).
"""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox

import tool_settings

MORE = "More: Help > Keyboard Shortcuts (F1) and docs/USER_GUIDE.md."

HELP = {
    "Sequence": [
        "Shows the bases, both strands, with features, primers and cut sites.",
        "Click to place the cursor; drag, or Shift+click, to select.",
        "Type A C G T N to edit; Ctrl+C / X / V copy, cut, paste (features come along).",
        "Click a feature or primer to select it; double-click to edit it.",
        "Right-click a feature, or selected bases, for primers, features and copy.",
        "Ctrl+wheel (or Ctrl+= / Ctrl+-) changes the letter size; Ctrl+0 resets it.",
        "Hover a feature, primer or cut for details.",
    ],
    "Map": [
        "Shows the whole sequence as a circle (plasmid) or a line, with its features.",
        "Click a feature to select it; double-click to edit it.",
        "Click or drag the backbone to place the cursor or select bases.",
        "Right-click a feature to edit, delete or show its translation.",
        "Hover a feature to see its name.",
        "Linear map: the wheel zooms at the pointer (or use the slider); Shift+wheel or middle-drag pans. Its zoom is separate from the Sequence view's.",
        "The selection is shared with the Sequence view.",
    ],
    "Alignment": [
        "Shows a stretch of the sequence like a genome browser: tracks, reads, primers, features.",
        "Wheel zooms at the pointer; drag pans; double-click zooms in.",
        "Shift+drag selects bases; click the overview bar to jump.",
        "Click a feature, primer or read to select it; right-click a track for its display.",
        "Tracks (header): show, hide, reorder, group and restyle many tracks at once; Fit squeezes them into the pane.",
        "Signal tracks (bigWig, bedGraph, wig): right-click > Signal for bars, a line or a heatmap, the range, overlays.",
        "Click a track name to hide the track; click a group's line to collapse it.",
        "Right-click a selection or the ruler: highlight it, recap every track there, export it as an image.",
        "Right-click a track's item: add it as a feature; a track: Copy to Features, Next Item; a signal: Regions Above Threshold.",
        "Type a position, range (12,000-15,000) or name in the locus box and press Enter.",
    ],
    "Features": [
        "Lists every feature: type, name, position, direction.",
        "Click a row to select the feature in all views.",
        "Delete removes the current feature.",
        "Right-click to delete, hide or add to the Features Database.",
    ],
    "Primers": [
        "Lists the primers from your database that bind this sequence, one group per primer.",
        "Click a row to select its binding site; double-click to edit the primer.",
        "New and Edit open the primer dialog; Database opens the full list.",
        "Right-click a primer for Edit, Hide in This File or Delete.",
        "The columns give length, Tm, mismatches and hairpin or dimer checks.",
        "Options sets the PCR conditions for Tm and the background genome.",
    ],
    "Protein": [
        "Shows the amino acids of the selected CDS, or of the selected DNA in the frame you choose.",
        "Click a CDS (or select DNA) in any view to show it.",
        "Length, mass, isoelectric point and charge are shown above the sequence.",
        "Select letters to select their codons in the DNA. The legend on the right explains and switches the tracks under them.",
        "Click a feature inside the CDS to select it.",
        "3D structure... gets a model (AlphaFold, PDB, ESMFold or your own command) and shows it beside the letters: drag to turn, click a residue to select it.",
        "Options chooses the genetic code (a CDS with /transl_table keeps its own).",
    ],
    "PCR": [
        "Left: the template, editable. Middle: the primers. Right: the product, read only.",
        "Pick a forward and a reverse primer that bind the template, or select a region on the left and press Design.",
        "The last choice (a primer in a list, or a designed pair) is the one shown on the right.",
        "Open Product makes it a new sequence; Save Pair to Primer Database keeps a designed pair.",
        "Click a primer on the template (then Shift+click one on the other strand) to choose it in the lists.",
    ],
    "Spacer Finder": [
        "Lists Cas9 / Cas12a guides (SpCas9, SaCas9 or Cas12a) with the PAM, on both strands.",
        "Choose the nuclease and the region, then Find guides.",
        "Click a row to select the site; the cut is drawn on the sequence.",
        "Check specificity counts other near-matches here and in the background genome.",
        "Add selected as feature annotates the guide.",
        "The background genome is set in the cog (Primers options).",
    ],
    "Local Align": [
        "Aligns a pasted sequence to the open one, both strands.",
        "Paste plain bases or FASTA, or use Paste from clipboard, then Align.",
        "Each matching piece is a row, in query order, with identity and what lies between pieces.",
        "Click a row to select its span in the sequence.",
        "Hits appear as a Query hits track in the Alignment view; Clear hits removes it.",
        "The Multiple alignment tab compares several sequences at once (the records of this file, files you pick, or the open sequence plus the selection): one row per sequence, columns coloured by conservation, a consensus row, and each sequence's identity to the reference or consensus.",
        "Click a column to put the cursor there; drag to select columns. A Terminal alignment (MAFFT, Clustal Omega) opens in this tab.",
    ],
    "Detect ORFs": [
        "Finds open reading frames on both strands, including across the origin of a circle.",
        "Set the minimum length at the top.",
        "Click a row to select the ORF in the sequence.",
        "Add selected as CDS features annotates the ones you pick.",
        "The list follows your edits.",
        "Options chooses the genetic code and whether GTG, TTG... can start an ORF.",
    ],
    "Digest": [
        "Shows where restriction enzymes cut and the fragments they give.",
        "Tick enzymes (the filter box searches names and sites).",
        "Click a fragment to select it in the sequence; hover a cut for details.",
        "The virtual gel shows each digest against a ladder.",
        "+ lane adds another digest to compare.",
    ],
    "Detect Features": [
        "Finds known parts (promoters, tags, resistance genes...) from the built-in and your own library.",
        "Tick the parts you want, then Add ticked features.",
        "Settings sets the minimum identity and length.",
    ],
    "Plasmid Checkup": [
        "Traces each transcription unit and says whether each CDS should be expressed.",
        "Verdicts: Expressed, Likely, Unlikely or Not expressed, with the reason.",
        "Select a row to select it in the sequence.",
        "Add as features puts the units on the map.",
        "It checks the design, not expression levels.",
        "Settings tunes the motif prediction.",
    ],
    "Assembly": [
        "Plans a Gibson, Golden Gate or Sidewinder assembly of up to 12 fragments.",
        "Click a fragment box to take its DNA from an open or saved file.",
        "Right-click a fragment for source, move, reverse complement and destination vector.",
        "Junction cards show whether the ends fit.",
        "Design primers adds the tails the method needs.",
        "Perform assembly opens the product as a new sequence.",
    ],
    "Recombineering": [
        "Replaces, deletes or inserts a region of a circular genome with a cassette.",
        "Genome: click its box, choose the file, select the region to replace, press Use.",
        "Cassette: choose its file, select the stretch to insert, press Use.",
        "Design primers gives the cassette the genome's homologies as 5' tails.",
        "Right-click the cassette to reverse complement it. Perform assembly opens the new genome.",
    ],
    "Synthesis": [
        "Will a synthesis company make this DNA quickly? Paste bases, or use the selection or the whole sequence.",
        "The Synthesis Success Calculator (Salis lab model) runs on this computer and needs ViennaRNA (Install ViennaRNA…).",
        "Each measure in its table says what it counts: long repeats, hairpins, GC-poor or GC-rich stretches, extreme melting temperatures.",
        "IDT's own complexity screen sends the sequence to idtdna.com, only when you press its button and confirm.",
    ],
}


def help_text(view):
    """The guide of `view` as plain text, one bullet per line, then where to find more."""
    return "\n".join(f"• {line}" for line in HELP.get(view, [])) + f"\n\n{MORE}"


def options_tip(view):
    """Tooltip of the cog: what it opens."""
    return f"Options of this view (Tools > Settings > {tool_settings.LABELS.get(view, view)})"


def open_options(window, view):
    """Open the Settings page of `view` (it says so when the view has no options)."""
    window.open_settings(view)


def show_view_help(parent, view):
    """The "?" popup: the view's guide and a Close button (Enter / Esc close it)."""
    box = QMessageBox(QMessageBox.NoIcon, f"{view}: Quick Guide", help_text(view), QMessageBox.Close, parent)
    box.setTextInteractionFlags(Qt.TextSelectableByMouse)
    box.exec()
