import discord
from random import randint
from gradio_client import Client

def stable_cascade(prompt):
    try:
        client = Client('multimodalart/stable-cascade')
        result = client.predict(
            prompt=prompt,
            negative_prompt='',
            seed=randint(0, 2147483647),
            width=1024,
            height=1024,
            prior_num_inference_steps=20,
            prior_guidance_scale=4,
            decoder_num_inference_steps=10,
            decoder_guidance_scale=0,
            num_images_per_prompt=1,
            api_name="/run",
        )
        return discord.File(result, 'image.webp')
    except ConnectionError as conn_err:
        print(f"Connection error: {conn_err}")
        return "Failed to connect to the server. Please try again later."
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return f"An error occurred: {e}"
