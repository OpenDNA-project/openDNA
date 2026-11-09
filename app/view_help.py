"""What the "?" and the cog next to a pane's view picker show, one entry per view of the picker (workspace.py: VIEWS and
TAB_TOOLS). Edit the texts here; tests/test_view_help.py checks that every view has one.

HELP[view]    - the short guide: 5-8 bullet lines, plain language (keep it near 120 words).
OPTIONS[view] - the page of Tools > Settings that holds the view's options (an empty page says it has none).
Every key gesture below is one the code handles (see also shortcuts.py, which lists them for Help > Keyboard Shortcuts).
"""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QLabel, QVBoxLayout

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
        "Click a track name to hide or show the track.",
        "Type a position, range (12,000-15,000) or name in the locus box and press Enter.",
    ],
    "Features": [
        "Lists every feature: type, name, position, direction.",
        "Click a row to select the feature in all views.",
        "Delete removes the current feature.",
        "Click the arrow beside a gene to show or hide its parts.",
        "Right-click to group as a gene, ungroup or add to the Features Database.",
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
        "Select letters to select their codons in the DNA. Colour colours them by property, hydrophobicity, size or charge.",
        "Click a feature inside the CDS to select it.",
        "3D structure... gets a model (AlphaFold, PDB, ESMFold or your own command) and shows it beside the letters: drag to turn, click a residue to select it.",
        "Options chooses the genetic code (a CDS with /transl_table keeps its own).",
    ],
    "PCR": [
        "Shows the product of two primers that bind this sequence.",
        "Pick a forward and a reverse primer from the lists.",
        "The product's size is shown at once.",
        "Open product makes it a new sequence.",
        "Tip: click one primer, then Shift+click one on the other strand in the Sequence view to select what they amplify.",
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
        "Parts found around a CDS can be grouped into a gene.",
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
}

# Where each view's options live: the page of the view's own name in Tools > Settings (its group is the picker's).
OPTIONS = {name: name for name in HELP}


def help_text(view):
    """The guide of `view` as plain text, one bullet per line, then where to find more."""
    return "\n".join(f"• {line}" for line in HELP.get(view, [])) + f"\n\n{MORE}"


def options_tip(view):
    """Tooltip of the cog: what it opens."""
    return f"Options of this view (Tools > Settings > {tool_settings.LABELS.get(view, view)})"


def open_options(window, view):
    """Open the Settings page of `view` (it says so when the view has no options)."""
    window.open_settings(view)
    return True


class ViewHelpDialog(QDialog):
    """The "?" popup: the view's guide and a Close button (Enter / Esc close it)."""

    def __init__(self, parent, view):
        super().__init__(parent)
        self.setWindowTitle(f"{view}: Quick Guide")
        self.label = QLabel(help_text(view), wordWrap=True)
        self.label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.label.setMinimumWidth(420)
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.Close).setDefault(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(buttons)
