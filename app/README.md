# CAT Rescoring Explorer (companion app)

Interactive Gradio app accompanying the paper. Lets readers explore the
simulation under user-specified parameters across four tabs:

1. Cross-population misalignment (ρ)
2. Random noise and structured shifts
3. Stepwise trajectory by test length
4. Single-point exploration of any ρ × test-length combination

**Live app:** https://huggingface.co/spaces/doriax/cat-rescoring-explorer

## Run locally

```bash
pip install gradio numpy matplotlib
python app.py
```

## Deploy on Hugging Face Spaces

Create a new Space with the **Gradio** SDK and add `app.py` plus a
`requirements.txt` containing:

```
gradio
numpy
matplotlib
```

The app re-implements the IRT/CAT simulation in NumPy (no R dependency), so it
runs standalone.
