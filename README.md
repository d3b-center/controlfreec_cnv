# KFDRC ControlFreeC WGS Workflow

## Summary
This workflow contains a wrapper around the Seven Bridges Genomics developed ControlFreeC v11.6 tool.  It creates the mini-pileup that
software would normally create internally of the b allele frequency file in parallel for improved performance.  It also renames outputs
and creates a bonus seg file output that is not native to ControlFreeC output

## Related Links/Info:
Tool Docker Pull: `images.sbgenomics.com/vojislav_varjacic/control-freec-11-6:v1`

Control FreeeC Docs: http://boevalab.inf.ethz.ch/FREEC/tutorial.html
